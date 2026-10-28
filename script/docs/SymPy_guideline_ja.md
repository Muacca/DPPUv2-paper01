# SymPy 実装ガイドライン — DPPUv2 エンジン

- **Target:** AI Coding Assistants & Contributors
- **Context:** 静的・一様な $S^3\times S^1$ 、 $T^3\times S^1$ 、コンパクト $\mathrm{Nil}^3\times S^1$ 上の、定数係数Nieh–Yan項を含むEinstein–Cartan作用

本ドキュメントは、このディレクトリの記号計算エンジン `DPPUv2_engine_core_v3.py` と三つの幾何runnerについて、実装上の指針を定める。添字順・フレーム・torsion・大域幾何の規約は[CONVENTIONS](CONVENTIONS_ja.md)で固定しており、論文では第2章と付録Aに同じ定義を示している。

⇒ [English version](SymPy_guideline.md) | [幾何規約](CONVENTIONS_ja.md)

-----

## 1\. 高速化エンジニアリング指針 (Optimization Rules)

SymPyを用いた積分計算における「処理時間爆発（2時間以上）」を防ぎ、数秒で完了させるための鉄則。

現在のrunnerは成分が一定の左不変フレームを用い、一様な密度に大域体積を掛けて作用を得る（ステップ `E4.12`）。以下の指針は、明示的な積分や大きな中間式を導入する場合に適用する。

### Rule 1.1: `expand()` + `cancel()` 戦略の徹底

積分（`integrate`）を行う直前に、被積分関数に対して高コストな `simplify()` を使用してはならない。代わりに「展開と約分」を行うこと。

  * **Don't:**
    ```python
    density = simplify(density)  # NG 禁止：巨大な式の因数分解は計算コストが極大
    result = integrate(density, x)
    ```
  * **Do:**
    ```python
    density = cancel(expand(density))  # OK 推奨：多項式の和にすることで項別積分を誘発
    result = integrate(density, x)
    ```

### Rule 1.2: 積分後の中間簡約の抑制

多重積分（例： $\phi$ 積分 → $\theta$ 積分）の際、中間結果に対して過度な `simplify` を行わない。`cancel` 程度に留め、最終的な積分直前で再度 `expand` する方が高速である。

-----

## 2\. テンソル成分の実装指針

### Rule 2.1: フレーム成分を所定の添字順で扱う

エンジンはEuclidean計量 $`g_{ab}=\delta_{ab}`$ （`metric_frame`、既定は `Matrix.eye(dim)`）の正規直交フレームで計算する。このときフレーム添字の上げ下げは恒等操作であり、成分を直接用いてよい。ただし配列の位置は[CONVENTIONS 第2節](CONVENTIONS_ja.md)の順序に従う。

* `C[a,b,c]` $`=C^a{}_{bc}`$
* `Gamma[a,b,c]` $`=\Gamma^a{}_{bc}`$ 。最後の添字が微分方向
* `Riemann[a,b,c,d]` $`=R^a{}_{bcd}`$

非正規直交基底や別の計量では、全ての添字を計量で下げてから入れ替え、必要に応じて上げ直す。配列の転置だけでは添字操作にならない。

### Rule 2.2: 入力配列 $\tau$ と幾何学的torsion $T$ を区別する

`torsion_tensor` には入力配列 $\tau$ を格納する。Cartan方程式から得られる幾何学的torsionは $T=-\tau$ であり、`geometric_torsion_tensor` に格納する。contortionのステップ `E4.5` はコードの式を $\tau$ に適用する。

```python
# Step E4.5: 入力配列tauのcontortion
for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            val = (tau[a, b, c] + tau[b, c, a] - tau[c, a, b]) * Rational(1, 2)
```

幾何学的torsionで書くと、同じcontortionは $`K_{abc}[T]=-\frac12(T_{abc}+T_{bca}-T_{cab})`$ である（論文の式(M2-03)）。配列と二つのパラメータを $`\tau=-T`$ 、 $`\eta_\tau=-\eta`$ 、 $`V_\tau=-V`$ と同時に変換する。 $K$ だけ、または一方のtorsionパラメータだけの符号を反転しない。

### Rule 2.3: 構造係数を含むtorsionの整合性検証

非座標フレームでは、torsion成分に構造係数が入る。

$$
T^a{}_{bc}=C^a{}_{bc}+\Gamma^a{}_{cb}-\Gamma^a{}_{bc}.
$$

EC接続を構築した後、これが $-\tau$ を再現し、ミスマッチが **0** であることを確認する。

```python
dim = engine.data['dim']
C = engine.data['structure_constants']
Gamma_EC = engine.data['connection_EC']
tau = engine.data['torsion_tensor']
mismatch = sum(
    1
    for a in range(dim) for b in range(dim) for c in range(dim)
    if simplify(C[a, b, c] + Gamma_EC[a, c, b] - Gamma_EC[a, b, c] + tau[a, b, c]) != 0
)
assert mismatch == 0
```

`Gamma_EC[a,b,c] - Gamma_EC[a,c,b]` だけでは $C$ が抜ける。 $C\ne0$ の $S^3$ やNilでは有効な検証にならない。

### Rule 2.4: エンジン内蔵の検証

* ステップ `E4.3` はLevi–Civita接続の計量適合性 $`\Gamma_{abc}+\Gamma_{bac}=0`$ を検証し、違反時は例外を送出する。
* ステップ `E4.7` は $`R_{abcd}=-R_{bacd}`$ と $`R_{abcd}=-R_{abdc}`$ を `PROVED_ZERO`・`WITNESS_NONZERO`・`UNPROVED` の判定で検証する。証明されていない残差を零として扱わない。
* torsionを持つ接続では、対交換 $`R_{abcd}=R_{cdab}`$ は要求しない（論文の式(B-02)）。

-----

## 3\. エンジンが用いる規約

### 3.1 Torsionと曲率

torsionと曲率は2-formで定義し、 $`\omega^a{}_b=\Gamma^a{}_{bc}e^c`$ とする（[CONVENTIONS 第6節](CONVENTIONS_ja.md)）。

$$
\begin{aligned}
T^a&=de^a+\omega^a{}_b\wedge e^b=\frac12T^a{}_{bc}\thinspace e^b\wedge e^c,\\
R^a{}_b&=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b=\frac12R^a{}_{bcd}\thinspace e^c\wedge e^d.
\end{aligned}
$$

### 3.2 接続

$$
\begin{aligned}
\Gamma^{\mathrm{EC}}_{abc}&=\Gamma^{\mathrm{LC}}_{abc}+K_{abc},\qquad
\Gamma^{\mathrm{LC}}_{abc}=\frac12\left(C_{abc}+C_{cba}-C_{bac}\right),\\
K_{abc}[T]&=-\frac12\left(T_{abc}+T_{bca}-T_{cab}\right)
=\frac12\left(\tau_{abc}+\tau_{bca}-\tau_{cab}\right).
\end{aligned}
$$

Levi–Civita部分は[CONVENTIONS 第5節](CONVENTIONS_ja.md)の一般Koszul公式であり、両側不変計量を仮定しない。

### 3.3 曲率成分

各runnerではフレーム成分が一定なので、ステップ `E4.7` は次の式を用いる。

$$
R^a{}_{bcd}=\Gamma^a{}_{ec}\Gamma^e{}_{bd}-\Gamma^a{}_{ed}\Gamma^e{}_{bc}+\Gamma^a{}_{be}C^e{}_{cd}.
$$

最後の項は非座標フレームの $de^c$ に由来し、省略してはならない（論文付録G）。

-----

## 4. Torsion Ansatz と Mode 分解規約

### 4.1 モード定義

| Mode | 入力成分 |
|---|---|
| `Mode.AX` | axial成分のみ（ $`V_\tau=0`$ ） |
| `Mode.VT` | vector-trace成分のみ（ $`\eta_\tau=0`$ ） |
| `Mode.MX` | 両成分 |

### 4.2 成分

ステップ `E4.4` は、runnerのパラメータ `eta` $`=\eta_\tau`$ 、`V` $`=V_\tau`$ とスケール `r` から入力配列を構成する。スケールの別名は $S^3$ で $r$ 、 $T^3$ で $`R_1`$ 、Nilで $R$ であり、ここでは $s$ と書く。

$$
\begin{aligned}
\tau^{\mathrm{ax}}_{ijk}&=\frac{2\eta_\tau}{s}\varepsilon_{ijk}\quad(i,j,k\in\lbrace 0,1,2\rbrace ),\\
\tau^{\mathrm{vec}}_{abc}&=\frac13\left(\delta_{ac}(V_\tau)_b-\delta_{ab}(V_\tau)_c\right),\qquad
(V_\tau)_a=(0,0,0,V_\tau).
\end{aligned}
$$

幾何学的torsion $T=-\tau$ は、 $`\eta=-\eta_\tau`$ 、 $`V=-V_\tau`$ として同じ形を持つ（論文の式(M2-05)）。runnerはステップ `E4.1`〜`E4.2` で計量・構造定数・体積・パラメータを与え、torsion配列を自前で構成しない。

-----

## 5. Nieh–Yan 密度と作用のバリアント

### 5.1 密度

Euclideanフレームで $`\varepsilon_{0123}=+1`$ とし、ステップ `E4.10` は次を計算する。

$$
N_{\mathrm{TT}}=\frac14\varepsilon^{abcd}T^e{}_{ab}T_{ecd},\qquad
N_{\mathrm{REE}}=\frac12\varepsilon^{abcd}R_{abcd},\qquad
N_{\mathrm{FULL}}=N_{\mathrm{TT}}-N_{\mathrm{REE}}.
$$

TTの係数 $1/4$ は二つのtorsion 2-formから、REEの係数 $1/2$ は一つの曲率2-formから生じる。 $`N_{\mathrm{TT}}`$ はtorsionの二次式なので、 $\tau$ で評価しても同じ値になる。三つのrunnerのansatzでは

$$
N_{\mathrm{TT}}=N_{\mathrm{REE}}=-\frac{4V\eta}{s}=-\frac{4V_\tau\eta_\tau}{s},\qquad N_{\mathrm{FULL}}=0
$$

である。

### 5.2 バリアント選択

| `NyVariant` | ラグランジアンのNY項 | 役割 |
|---|---|---|
| `NyVariant.FULL` | $`\theta N_{\mathrm{FULL}}`$ | 完全NY作用 |
| `NyVariant.TT` | $`\theta N_{\mathrm{TT}}`$ | 診断用作用 |
| `NyVariant.REE` | $`+\theta N_{\mathrm{REE}}`$ | 正符号で加える診断用作用。FULLの負の成分ではない |

このansatzではTTとREEが同じポテンシャルを与える。それらの障壁や $\theta$ 依存性は、完全NY作用の結果ではない。

### 5.3 実装

`step_E4_10_nieh_yan_frame` は `ny_density_TT`・`ny_density_Ree`・`ny_density_full` を格納し、`BaseFrameEngine` の引数 `ny_variant` で採用する密度を選ぶ。

-----

## 6. 作用と縮約ポテンシャル

### 6.1 作用

$$
L=\frac{\mathcal R_{\mathrm{EC}}}{2\kappa^2}+\theta N,\qquad
S=L\times\mathrm{Vol},\qquad
U=-S.
$$

ステップ `E4.11`・`E4.12`・`E4.13` で $L$ 、 $S$ 、 $U$ を構成する。体積 `total_volume` は各runnerが定める（[CONVENTIONS 第9節](CONVENTIONS_ja.md)）。 $U$ は `potential` に格納され、ログでは「Effective Potential V(r)」と表示される。

### 6.2 ポテンシャルの取得

```python
from DPPUv2_engine_core_v3 import Mode, NyVariant
from DPPUv2_runner_S3S1_v3 import S3S1Engine

engine = S3S1Engine(Mode.MX, NyVariant.FULL)
engine.run()
U = engine.get_effective_potential_function()
# 引数: (scale, V_tau, eta_tau, theta_NY, L, kappa)
```

* スケール引数は $S^3$ で `r`、 $T^3$ で `R1`、Nilで `R` である。
* torsionの引数はコード入力である。論文の幾何学的な値を使う場合は $`V_\tau=-V`$ 、 $`\eta_\tau=-\eta`$ を渡す。
* $T^3$ では保存されたポテンシャルに `R2`・`R3` が残る。`get_effective_potential_function()` を呼ぶ前に、変分経路（例： $`R_2=\alpha R_1`$ 、 $`R_3=\beta R_1`$ ）を `engine.data['potential']` に代入する。
* `get_potential_decomposition()` はスケールの冪ごとの係数を返す。

-----

## 7. 半径極値の分類

ステップ `E4.14` は微分を記録するだけである。分類には `DPPUv2_parameter_scan_v3.py` の解析的分類器を用いる。

* `get_cubic_coefficients(topology, ny_variant, V, eta, theta)` は、登録済みの球面・等方 $T^3$ ・一様Nilの経路について、コード入力で $U=ar^3+br^2+cr$ の $(a,b,c)$ を返す。
* `analyze_cubic_stability(a, b, c)` は正半径領域で分類する（論文第4章・表3）。

| 分類 | $r>0$ での定義 |
|---|---|
| Type-I | 正半径の局所極大の後に厳密な局所極小がある |
| Type-II | 先行する正半径の極大を伴わない厳密な局所極小がある |
| Type-III | 厳密な局所極小がない（単調・平坦・runaway・停留変曲点） |

* type-Iの高さは極大から極小までの作用差、type-IIの深さは原点極限 $U(0)=0$ からの差である。同じ種類の障壁として描画しない。
* $`b^2-3ac\le64\epsilon_{\mathrm{float}}\max(b^2,\lvert 3ac\rvert)`$ の場合は、厳密な極小を解像しないとして扱う。判別式が零であることの証明ではなく、解析的不等式を優先する。
* 三次係数を取得できないcallableは有界な数値探索へfallbackし、その結論は探索窓と分解能に制限される。

これらは固定した変分経路上のoff-shellな半径極値である。全Einstein–Cartan真空、動的安定性、核生成、物理的な質量や崩壊率を示すものではない（論文第5章・第7章）。

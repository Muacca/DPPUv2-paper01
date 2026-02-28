# SymPy 実装ガイドライン — dppu エンジン

- **Target:** AI Coding Assistants & Contributors
- **Context:** Einstein-Cartan Theory on Curved Spacetime (S³×S¹, T³×S¹, Nil³×S¹)

本ドキュメントは、`dppu` パッケージにおける数式処理（SymPy）の実装指針、および理論物理学的規約（Convention）を定めるものである。以下のルールを厳守すること。

⇒ [English version](SymPy_guideline.md) | [幾何規約](CONVENTIONS_ja.md)

-----

## 1\. 高速化エンジニアリング指針 (Optimization Rules)

SymPyを用いた積分計算における「処理時間爆発（2時間以上）」を防ぎ、数秒で完了させるための鉄則。

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

## 2\. 理論実装指針 (Theoretical Implementation Rules)

曲がった時空（Curved Spacetime）において、テンソル演算の整合性を保つための鉄則。

### Rule 2.1: 配列インデックス操作の禁止 (Robust Method)

非対角計量や $g_{\mu\nu} \neq 1$ の環境下では、`T[mu, nu, lam]` のように配列のインデックスを入れ替える操作は、物理的なテンソルの添字操作（上げ下げ）と等価ではない。

  * **Don't:**
    ```python
    # NG 禁止：物理的に誤った成分を参照する恐れがある
    term = T_tensor[mu, nu, lam]
    ```
  * **Do:** 必ず計量 $g_{\mu\nu}$ を介して操作する。
    1.  全ての添字を下げて完全共変形 $T_{\lambda\mu\nu}$ を作る。
    2.  添字の入れ替え（Permutation）を行う。
    3.  必要に応じて計量で添字を上げる。

#### 正規直交フレーム基底での最適化

正規直交フレーム基底においては、計量が単位行列（ $g_{ab} = \eta_{ab} = \text{diag}(1,1,1,1)$ または $\text{diag}(-1,1,1,1)$）であるため、添字の上げ下げ計算を省略し、直接成分計算を行うことで高速化する。

ただし、物理的定義（Hehl 1976）の符号パターン $(+1, +1, -1)$ は厳守すること。

```python
# ============================================================
# Golden Logic for Contortion (Frame Basis / DPPUv2 Standard)
# ============================================================
# Assumption: Metric is diagonal/identity (Orthonormal Frame)
# Therefore T^a_bc and T_abc behave identically in code logic.

K_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            # Formula: K_abc = (1/2)(T_abc + T_bca - T_cab)
            # Note: Using T[a,b,c] directly as T_abc
            
            term = (T_tensor[a, b, c] + T_tensor[b, c, a] - T_tensor[c, a, b])
            
            val = term * Rational(1, 2)
            
            if val != 0:
                K_tensor[a, b, c] = cancel(expand(val))
```

### Rule 2.2: 自己無撞着性の検証 (Consistency Check)

EC接続 ($\Gamma_{\text{EC}}$) を構築した後、必ず以下の検証コードを実行し、ミスマッチが **0** であることを確認しなければならない。

```python
# Torsion Consistency Check
T_verify = Gamma_EC[lam, mu, nu] - Gamma_EC[lam, nu, mu] # Hehl定義
mismatch = count(simplify(T_verify - T_original) != 0)
assert mismatch == 0
```

-----

## 3\. 理論的規約 (Standard Conventions)

論文執筆時の混乱を防ぐため、**Hehl (1976) 標準** に準拠する。

### 3.1 Torsion Definition

捩れテンソル $T^\lambda_{\ \mu\nu}$ の定義：
$$T^\lambda_{\ \mu\nu} \equiv \Gamma^\lambda_{\ \mu\nu} - \Gamma^\lambda_{\ \nu\mu}$$
（注：本エンジンのフレーム基底での捩れ成分は、CONVENTIONS の第6節で定義された捩れ2-form $T^a = de^a + \omega^{a}{}\_b\wedge e^b$ から $T^a = \frac{1}{2}T^{a}{}\_{bc}\,e^b\wedge e^c$ の係数比較により抽出される）

### 3.2 Contortion Formula

上記のTorsion定義と整合するContortion $K^\lambda_{\ \mu\nu}$ の公式（Verified Formula）：

$$K_{\lambda\mu\nu} = \frac{1}{2} \left( T_{\lambda\mu\nu} + T_{\mu\nu\lambda} - T_{\nu\lambda\mu} \right)$$

  * 符号パターン: **$(+1, +1, -1)$**
  * 注意: この式は $T_{\lambda\mu\nu}$（全ての添字を下げたもの）に対して適用すること。

### 3.3 Einstein-Cartan Connection

$${\Gamma_{\text{EC}}}^\lambda_{\ \mu\nu} = {\Gamma_{\text{LC}}}^\lambda_{\ \mu\nu} + K^\lambda_{\ \mu\nu}$$

  * $\Gamma_{\text{LC}}$: Levi-Civita接続（Christoffel記号）
  * $K$: Contortion

-----

## 4. Torsion Ansatz と Mode 分解規約

$M^3 \times S^1$ ミニスーパースペース・アンサーツでの捩れテンソルは、以下の3モードで指定する。

### 4.1 モード定義

| Mode | 物理的成分 | パラメータ |
|---|---|---|
| `Mode.AX` | 軸対称成分（T1）のみ | $\eta \neq 0$, $V = 0$ |
| `Mode.VT` | ベクトル跡部分（T2）のみ | $\eta = 0$, $V \neq 0$ |
| `Mode.MX` | T1 + T2 の両成分 | $\eta \neq 0$, $V \neq 0$ |

### 4.2 物理的対応

- **T1（軸対称成分）**: 軸性ベクトル $S^\mu = (\eta/r)(0,0,0,1)$ の双対。空間添字 $a,b,c \in \{0,1,2\}$ に対して $T_{abc} = (2\eta/r)\,\varepsilon_{abc}$。
- **T2（ベクトル跡部分）**: ベクトル $V_\mu = V\,\delta^3_\mu$（$\tau$ 成分のみ）の双対。 $T_{abc} = \frac{1}{3}(\delta_{ac}V_b - \delta_{ab}V_c)$。

### 4.3 実装ルール

`dppu/torsion/ansatz.py` の `construct_torsion_tensor(mode, r, eta, V, metric, dim)` を使用すること。
$T_{abc}$ を手打ちで構築してはならない。

-----

## 5. Nieh-Yan トポロジカル項のバリアント

### 5.1 Nieh-Yan 分解

完全な Nieh-Yan 密度：
$$N = N_{\mathrm{TT}} - N_{\mathrm{Ree}},$$
$$N_{\mathrm{TT}} = \frac{1}{4}\varepsilon^{abcd}T^{e}{}\_{ab}T_{ecd},\qquad N_{\mathrm{Ree}} = \frac{1}{4}\varepsilon^{abcd}R_{abcd}.$$

### 5.2 バリアント選択

| `NyVariant` | 使用する密度 |
|---|---|
| `NyVariant.TT` | $N_{\mathrm{TT}}$ のみ |
| `NyVariant.REE` | $N_{\mathrm{Ree}}$ のみ |
| `NyVariant.FULL` | $N_{\mathrm{TT}} - N_{\mathrm{Ree}}$（標準） |

### 5.3 実装

パイプラインステップ `E4.10` で全バリアントが計算される。バリアント選択は engine の `__init__` で `ny_variant` 引数を指定する。`dppu/torsion/nieh_yan.py` を参照。

-----

## 6. 拡張ラグランジアンと Weyl 結合定数 $\alpha$

### 6.1 作用の形式

$$S = \int L \times \mathrm{Vol},\qquad
L = \frac{R}{2\kappa^2} + \theta_{\mathrm{NY}}\times N + \alpha\times C^2.$$

| パラメータ | 意味 |
|---|---|
| $\kappa$ | Einstein-Cartan 重力結合定数 |
| $\theta_{\mathrm{NY}}$ | Nieh-Yan 結合定数（トポロジカル） |
| $\alpha$ | Weyl 結合定数（共形不変項） |

$\alpha \leq 0$ では定理1により安定真空が保護される。 $\alpha > 0$ では定理2により有効ポテンシャルが非有界となる。

### 6.2 有効ポテンシャルの取得

```python
engine.run()
V_func = engine.get_effective_potential_function()
# 呼び出しシグネチャ：
# V_func(r, V_param, eta, theta_NY, L, kappa, epsilon, alpha) -> float
```

$V_{\rm eff} = -S$。パイプラインステップ `E4.13` で抽出される。

### 6.3 実装位置

`dppu/action/lagrangian.py` の `compute_lagrangian()`、`dppu/action/potential.py` を参照。

-----

## 7. 数値最適化戦略（Phase Atlas 探索）

### 7.1 2段階戦略

**Stage 1：Brute-force グリッド探索**

`scipy.optimize.brute`（`Ns` 点/軸）で $(r, \varepsilon)$ 2D グリッドを粗く探索し、大域極小の領域を特定する。

**Stage 2：Multi-start L-BFGS-B 精密化**

Stage 1 の上位 $N$ 候補を初期点として `scipy.optimize.minimize`（L-BFGS-B、`ftol=1e-8`）を実行し、高精度最小値を得る。

### 7.2 安定性分類

$(r^\*, \varepsilon^\*)$ の最小値発見後、以下に分類する：

| 分類 | 条件 | 物理的意味 |
|---|---|---|
| Type-I | $V(r)$ が $r \to 0$ で立ち上がる（障壁あり） | 安定真空（核生成障壁あり） |
| Type-II | $V(r)$ が $r \to 0$ で単調減少（障壁なし） | 自発的核生成が可能 |
| Type-III | 物理的領域に局所最小なし | 不安定配位 |

`dppu/action/stability.py` の `analyze_stability()` を使用すること。

### 7.3 注意事項

- $\alpha > 0$ の場合、最適化が探索境界 $(r \to 0^+,\,\varepsilon \to -1^+)$ に張り付く（`converged = False`）。これは最適化の失敗ではなく、ポテンシャルが探索範囲内で非有界であることの正確な報告である。
- $\mathrm{Nil}^3$ の平坦極限確認のため、探索範囲を通常より広く $\varepsilon_{\rm max} = 5.0$ に設定している（`scripts/paper02/` 参照）。

-----

## 8. Self-duality（SD）診断規約

### 8.1 曲率の Hodge 双対

$$(*R)^{ab}{}_{cd} = \frac{1}{2}\varepsilon_{cdef}\,R^{ab,ef},$$

ここで $\varepsilon_{cdef}$ はフレーム基底の Levi-Civita 記号（`dppu/utils/levi_civita.py` の `levi_civita_4d()` を使用）。

### 8.2 Pontryagin 内積と SD 残差

$$E_{RR} = \langle R, R\rangle = R_{abcd}R^{abcd},\qquad
P = \langle R, *R\rangle = R_{abcd}(*R)^{abcd}.$$

| 条件 | 物理的状態 |
|---|---|
| SD 残差 $< \varepsilon_{\rm SD}$ かつ $\|R\| > \varepsilon_R$ | Self-dual |
| ASD 残差 $< \varepsilon_{\rm SD}$ かつ $\|R\| > \varepsilon_R$ | Anti-self-dual |
| $P = 0$（命題1） | Chiral equilibrium（ $M^3 \times S^1$ 上で代数的に成立） |

### 8.3 使用方法

```python
from dppu.curvature.self_duality import SDExtensionMixin
SDExtensionMixin.attach_to(engine)          # メソッドを動的にアタッチ
R = engine.get_R_ab_cd_numerical(params_dict)
diag = engine.evaluate_sd_status(params_dict)
# diag['P_RstarR'] で Pontryagin 値を確認（== 0 が期待値）
```

`dppu/curvature/self_duality.py`、`dppu/curvature/pontryagin.py` を参照。



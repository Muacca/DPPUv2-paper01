# CONVENTIONS — dppu エンジン幾何規約

本書は、`BaseFrameEngine` と各 topology runner が共有する **幾何・添字・符号規約**を固定する。
**全 runner はこの規約に従って `metric_frame` と `structure_constants` を定義すること。**

⇒ [English version](CONVENTIONS.md) | [SymPy ガイドライン](SymPy_guideline_ja.md)

## 1. 作用域と前提

* ここで扱う量はすべて **フレーム（正規直交基底）** 上の成分で表す。
* `metric_frame` はフレーム計量 $(g_{ab})$ であり、既定では $(g_{ab}=\delta_{ab})$（`Matrix.eye(dim)`）。
* DPPUv2 engine の曲率成分計算は、（現状）**フレーム方向微分が不要になる状況**を前提にしている。
  具体的には、構造定数 $(C^{a}{}\_{bc})$ と接続係数 $(\Gamma^{a}{}\_{bc})$ が「フレームに関して定数扱い」になる（左不変フレームなど）設定を runner が採用すること。
  * **Note:** この設計は $S^3 \times S^1$（Lie群）や Nil 多様体のような等質空間を扱う上では**合理的かつ効率的**である。ただし、将来的に対称性の低い一般的な曲率時空を扱う場合は、この前提がボトルネックになり得る点に留意すること。

## 2. 添字・配列のインデックス順

配列格納は以下で固定する：

* 構造定数：`C[a,b,c] = C^a_{bc}`
* 接続係数：`Gamma[a,b,c] = Γ^a_{bc}`
* リーマン曲率：`Riemann[a,b,c,d] = R^a_{bcd}`

添字の意味：

* $(a)$：出力（上付き）成分
* $(b,c,d)$：入力（下付き）成分
  とくに $(\Gamma^{a}{}\_{bc})$ は $(\nabla_{E_c} E_b = \Gamma^{a}{}\_{bc} E_a)$ に対応する（最後の $(c)$ が “微分方向”）。

## 3. フレーム・コフレームと構造定数の定義（ここが最重要）

フレームの双対を $(\{E_a\})$、コフレーム（1-forms）を $(\{e^a\})$ とする。

### 3.1 コフレームの構造方程式（固定）

$$
de^a = \frac12 C^a{}_{bc} e^b\wedge e^c,
\qquad C^a{}_{bc} = - C^a{}_{cb}.
$$

### 3.2 双対フレームの交換関係（同値）

上の定義は次と同値：

$$
[E_b, E_c] = - C^{a}{}\_{bc} E_a.
$$

> 注意：多くの教科書では $([E_b,E_c]=+f^{a}{}\_{bc}E_a)$ を採用する。
> 本プロジェクトでは **その $(f^{a}{}\_{bc})$ に対して $(C^{a}{}\_{bc}=-f^{a}{}\_{bc})$** の規約を採用している。

### 3.3 runner 実装ルール（推奨）

* **C は手打ちしない**。可能なら runner 側で $de^a$ を明示し、係数比較で $C^a_{bc}$ を抽出して `self.data['structure_constants']=C` に入れる。
* 最低限、 $C^{a}_{bc}$ が **b,c で反対称**になっていることを自動チェックすること。

## 4. 接続（スピン接続）とメトリック整合

接続 1-form を

$$
\omega^{a}{}\_{b} = \Gamma^{a}{}\_{bc} e^c
$$

で定義する。

メトリック整合（ローレンツ接続／直交接続）を仕様として固定：

$$
\omega_{ab} = -\omega_{ba}
\quad(\Leftrightarrow\quad
\Gamma_{abc} = -\Gamma_{bac})
$$

ただし $(\Gamma_{abc} = g_{ad}\Gamma^{d}{}\_{bc})$。

## 5. Levi-Civita 接続（DPPUv2 engine の一般 Koszul 実装）

フレームが正規直交で、上記の構造定数規約を採用したとき、Levi-Civita 接続は DPPUv2 engine では次の**一般 Koszul 公式**で計算する：

$$
\Gamma^a{}_{bc}
= \frac12\Big(
C^a{}_{bc} + C^c{}_{ba} - C^b{}_{ac}
\Big).
$$

（これは本書 3.2 の交換関係の符号を採用した場合の形である。）

**重要な注記:**

1. この公式は **bi-invariant 計量を仮定しない**。
   左不変フレーム上の Levi-Civita 接続として、Nil³ のような非 bi-invariant な場合でも正しく機能する。

2. SU(2) のように（低い添字で）構造定数が全反対称 $C_{abc} = -C_{bac} = -C_{acb}$ になる特殊な場合は、
   この公式は $\Gamma^a_{bc} = \frac{1}{2} C^a_{bc}$ に簡約される。

3. engine は計算後に **metric compatibility** $\Gamma_{abc} + \Gamma_{bac} = 0$ を自動検証する。
   これに違反する場合は実装エラーとして即座に例外を投げる。

## 6. ねじれと曲率

ねじれ 2-form：

$$
T^a = de^a + \omega^{a}{}\_b \wedge e^b,
\qquad
T^a = \frac12 T^{a}{}\_{bc} e^b\wedge e^c.
$$

曲率 2-form：

$$
R^a{}_b = d\omega^{a}{}\_b + \omega^{a}{}\_c\wedge \omega^{c}{}\_b,
\qquad
R^a{}_b = \frac12 R^{a}{}\_{bcd} e^c\wedge e^d.
$$

## 7. 曲率成分の計算式（engine が実際に使う形）

DPPUv2 engine の $R^{a}_{bcd}$ は（現状）次の形を用いる：

$$
R^{a}{}\_{bcd} = \Gamma^{a}{}\_{ec}\Gamma^{e}{}\_{bd} -\Gamma^{a}{}\_{ed}\Gamma^{e}{}\_{bc} +\Gamma^{a}{}\_{be} C^{e}{}\_{cd}.
$$

> 重要：一般にはここにフレーム方向微分項
> $(E_c(\Gamma^{a}{}\_{bd}) - E_d(\Gamma^{a}{}\_{bc}))$
> が現れるが、DPPUv2 engine ではそれを明示的に扱っていない。
> したがって runner は、左不変フレーム等により **$(\Gamma)$ がフレーム方向で定数扱い**になる設定を採用すること。

## 8. 必須セルフチェック（runner が満たすべき整合性）

runner は以下を満たすこと（落ちたら定義が engine と不整合）：

1. 構造定数の反対称：

$$
C^{a}{}\_{bc} + C^{a}{}\_{cb} = 0.
$$

2. メトリック整合（直交接続）：

$$
\omega_{ab} + \omega_{ba} = 0.
$$

3. リーマンの反対称（engine の strict check 対象）：

$$
R_{ab cd} = -R_{ba cd},\qquad
R_{ab cd} = -R_{ab dc}.
$$

---

## 9. Weyl テンソルと共形スカラー

### 9.1 Weyl テンソルの定義（4次元）

$$
C_{abcd} = R_{abcd} - \frac{1}{2}(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac}) + \frac{R}{6}(g_{ac}g_{bd} - g_{ad}g_{bc}).
$$

主要な性質：
- **無跡**： $C^{a}{}\_{bad} = 0$（全添字対で成立）
- **共形不変**：計量の共形変換 $g_{ab} \to \Omega^2 g_{ab}$ のもとで不変
- **共形平坦判定**： $C_{abcd} = 0 \Leftrightarrow$ 共形的平坦

フレーム基底（正規直交）での注意： $g_{ab} = \delta_{ab}$ のため添字の上げ下げは恒等変換に等しく、
 $C_{abcd} = C^{abcd}$ として直接成分を扱える。

### 9.2 Weyl スカラー

$$
C^2 = C_{abcd}\,C^{abcd} = \sum_{a,b,c,d} C_{abcd}^2.
$$

フレーム基底での高速計算：添字上げ下げを省略し、直接2乗和をとる。
$C^2 = 0 \Leftrightarrow$ 共形的平坦（等方的 $S^3 \times S^1$ で $\varepsilon = 0$ のとき成立）。

### 9.3 engine での実装位置

- モジュール：`dppu/curvature/weyl.py`
  - `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)` → $C_{abcd}$
  - `compute_weyl_scalar(C_abcd, metric_inv, dim)` → $C^2$
- パイプラインステップ：`E4.3b`（Levi-Civita 曲率計算の直後）

---

## 10. Squashed 等質空間と $\varepsilon$-パラメータ

### 10.1 Squashing の定義

体積保存の異方性変形パラメータ $\varepsilon$ を導入する。
$\varepsilon = 0$ が等方（アイソトロピック）基準点。物理的範囲： $\varepsilon \in (-1, +\infty)$（ $\varepsilon = -1$ で構造定数が発散し特異点）。

### 10.2 トポロジー別の構造定数スケーリング

左不変フレームの基底構造定数 $C^a{}_{bc}(\varepsilon=0)$ を以下のスケーリングで変形する。

**$S^3 \times S^1$（SU(2)）:**

| フレーム添字 $a$ | スケール因子 $\lambda_a(\varepsilon)$ |
|---|---|
| $a \in \{0, 1\}$ | $(1+\varepsilon)^{2/3}$ |
| $a = 2$ | $(1+\varepsilon)^{-4/3}$ |

体積保存の確認： $\lambda_0\lambda_1\lambda_2 = (1+\varepsilon)^{2/3+2/3-4/3} = 1$。

**$\mathrm{Nil}^3 \times S^1$（Heisenberg 群）:**

| 非自明な添字 $a$ | スケール因子 |
|---|---|
| $a = 2$（非可換成分） | $(1+\varepsilon)^{-4/3}$ |

$\varepsilon \to +\infty$ で構造定数が消失し、平坦 $T^3$ 的挙動に漸近する。

**$T^3 \times S^1$（Abelian 群）:**

構造定数は恒等的にゼロ（ $C^a{}_{bc} = 0$）。 $\varepsilon$ による変形は定義されない。
$C^2 = 0$ が全域で成立し、Weyl 項は常にゼロとなる（null test）。

### 10.3 物理的極限

| 極限 | 物理的意味 |
|---|---|
| $\varepsilon = 0$ | 等方的 $S^3$ ： $C^2 = 0$、Paper I の安定真空 |
| $\varepsilon \to +\infty$ | $\mathrm{Nil}^3$ の平坦極限（ $C^2 \to 0$ 漸近） |
| $\varepsilon \to -1^+$ | $S^3$ 構造の特異点（物理的到達不能） |
| $\varepsilon = -2$ | $C^2 = 0$ の数学的根だが $\varepsilon < -1$ のため物理的除外 |


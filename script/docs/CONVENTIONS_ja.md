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

### 6.1 入力配列と幾何学的torsion

`torsion_tensor` は入力パラメータ配列τを格納する。幾何学的成分は
`T_geom[a,b,c] = C[a,b,c] + Gamma[a,c,b] - Gamma[a,b,c] = -tau[a,b,c]`。
`geometric_torsion_tensor` はこの成分を公開する。contortionは
`K[a,b,c]=(tau[a,b,c]+tau[b,c,a]-tau[c,a,b])/2`で与える。

### 6.2 Nieh–Yan密度

Euclideanフレーム計量とepsilon(0,1,2,3)=+1の下で、成分縮約は
`N_TT=epsilon*T*T/4`、`N_REE=epsilon*R/2`、`N_FULL=N_TT-N_REE`である。
REEは作用に`+theta*N_REE`を加える診断量を表す。
3種のrunnerの静的一様ansatzでは`TT=REE=-4*V*eta/r`（各runnerのscale alias）、`FULL=0`。

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
| $\varepsilon = 0$ | 等方的 $S^3$ の背景： $C^2 = 0$ |
| $\varepsilon \to +\infty$ | $\mathrm{Nil}^3$ の平坦極限（ $C^2 \to 0$ 漸近） |
| $\varepsilon \to -1^+$ | $S^3$ 構造の特異点（物理的到達不能） |
| $\varepsilon = -2$ | $C^2 = 0$ の数学的根だが $\varepsilon < -1$ のため物理的除外 |

## 11. 大域幾何と積分領域

3種類のtopology runnerで使う大域背景を定める。計量は正規直交coframeの平方和で、スケールは全て正とする。共通の円について、$u\in[0,1)$、$e^3=L\,du$とし、**Lは円周長**を表す。engineは一様なスカラー密度に、指定した領域のRiemann体積を掛けて積分する。

### 11.1 丸いS³

$C^i{}_{jk}=4\varepsilon_{ijk}/r$に現れるスケールrは丸球半径の2倍であり、実半径は$a=r/2$。明示的なSU(2) coframeは$e^i=-r\lambda^i/4$で、

$$
\begin{aligned}
\lambda^0&=\cos\psi\,d\vartheta+\sin\psi\sin\vartheta\,d\varphi,\\
\lambda^1&=-\sin\psi\,d\vartheta+\cos\psi\sin\vartheta\,d\varphi,\\
\lambda^2&=d\psi+\cos\vartheta\,d\varphi.
\end{aligned}
$$

積分域は$0\le\vartheta\le\pi$、$0\le\varphi<2\pi$、$0\le\psi<4\pi$とし、極での退化と境界の同定$(\varphi+2\pi,\psi)\sim(\varphi,\psi+2\pi)$を伴う。3方向独立周期のトーラス座標ではない。正の計量密度は$r^3\sin\vartheta/64$。上のcoframeの向きは$d\vartheta\wedge d\varphi\wedge d\psi$と逆になる。

$$
\operatorname{Vol}(S^3)=\frac{\pi^2r^3}{4}=2\pi^2a^3,\qquad
\operatorname{Vol}(S^3\times S^1)=\frac{\pi^2Lr^3}{4},\qquad
R_{\mathrm{LC}}=\frac{24}{r^2}=\frac{6}{a^2}.
$$

### 11.2 T³

$i=1,2,3$について、$x^i\in[0,2\pi)$、$e^{i-1}=R_i\,dx^i$とする。$R_i$はcoframeのスケール係数で、空間各方向の周長は$\ell_i=2\pi R_i$。

$$
\operatorname{Vol}(T^3\times S^1)=(2\pi)^3LR_1R_2R_3
=L\ell_1\ell_2\ell_3.
$$

scanは$R_1=R_2=R_3=r$を使う。固定形状経路$R_2=\alpha r$、$R_3=\beta r$では正の$\alpha,\beta$を固定する。R2・R3自体を固定する場合は別の変分経路となる。

### 11.3 コンパクトNil³

Heisenberg群の積と左不変coframeを

$$
\begin{aligned}
(x,y,z)(x',y',z')&=(x+x',y+y',z+z'+xy'),\\
e^0&=R\,dx,\qquad e^1=R\,dy,\qquad e^2=R(dz-x\,dy)
\end{aligned}
$$

とする。コンパクト空間には左商$\Gamma\backslash H$を使い、格子を

$$
\Gamma=\{(2\pi j,2\pi k,4\pi^2l):j,k,l\in\mathbb Z\}
$$

と定める。生成元による座標同定は

$$
\begin{aligned}
(x,y,z)&\sim(x+2\pi,y,z+2\pi y),\\
(x,y,z)&\sim(x,y+2\pi,z),\\
(x,y,z)&\sim(x,y,z+4\pi^2).
\end{aligned}
$$

各同定はcoframeを保つ。周期$(a,b,c)=(2\pi,2\pi,4\pi^2)$に対し、水平生成元の交換子は中央方向の$ab=c$だけの並進であり、$ab/c=1$となる。基本領域は$[0,2\pi)\times[0,2\pi)\times[0,4\pi^2)$で、x面は上記のねじれを伴って貼り合わせる。正の密度$R^3\,dx\,dy\,dz$で積分する。

$$
de^2=-\frac1R e^0\wedge e^1,\qquad C^2{}_{01}=-\frac1R,\qquad
\operatorname{Vol}(\Gamma\backslash H\times S^1)=(2\pi)^4LR^3.
$$

底のトーラスの周長は各$2\pi R$、中央S¹繊維の周長は$4\pi^2R$。Rを変分するときも格子とその形状を固定する。この商は局所構造定数と指定体積を同時に実現し、中央座標の周期は$2\pi$ではない。

### 11.4 入力パラメータと実長さ

runner・消費側のAPIは局所スケールr/R1/R2/R3/Rと入力配列のパラメータ$\eta_\tau,V_\tau$を使用する。幾何学的torsionとの辞書は$\eta_g=-\eta_\tau$、$V_g=-V_\tau$。

同じS³のaxial振幅を$2\eta_g/r$から$2\chi/a$で表す場合は、$a=r/2$、$\chi=\eta_g/2$を使う。T³を周長で表す場合は、$2\zeta/\ell_1=2\eta_g/R_1$より$\zeta=2\pi\eta_g$とする。入力区間、刻み、極値もこの辞書に従って変換する。Vは符号付き係数のままで、追加の倍率を付けない。

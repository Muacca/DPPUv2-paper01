# CONVENTIONS (DPPUv2 v3)

本書は、`BaseFrameEngine (v3)` と各 topology runner が共有する **幾何・添字・符号規約**を固定する。
**全 runner はこの規約に従って `metric_frame` と `structure_constants` を定義すること。**

⇒ [English version](DPPUv2_engine_core_v3_CONVENTIONS.md)

## 1. 作用域と前提

* ここで扱う量はすべて **フレーム（正規直交基底）** 上の成分で表す。
* `metric_frame` はフレーム計量 $(g_{ab})$ であり、既定では $(g_{ab}=\delta_{ab})$（`Matrix.eye(dim)`）。
* v3 engine の曲率成分計算は、（現状）**フレーム方向微分が不要になる状況**を前提にしている。
  具体的には、構造定数 $(C^a{}_{bc})$ と接続係数 $(\Gamma^a{}_{bc})$ が「フレームに関して定数扱い」になる（左不変フレームなど）設定を runner が採用すること。
  * **Note:** この設計は $S^3 \times S^1$（Lie群）や Nil 多様体のような等質空間を扱う上では**合理的かつ効率的**である。ただし、将来的に対称性の低い一般的な曲率時空を扱う場合は、この前提がボトルネックになり得る点に留意すること。

## 2. 添字・配列のインデックス順

配列格納は以下で固定する：

* 構造定数：`C[a,b,c] = C^a_{bc}`
* 接続係数：`Gamma[a,b,c] = Γ^a_{bc}`
* リーマン曲率：`Riemann[a,b,c,d] = R^a_{bcd}`

添字の意味：

* $(a)$：出力（上付き）成分
* $(b,c,d)$：入力（下付き）成分
  とくに $(\Gamma^a{}_{bc})$ は $(\nabla_{E_c} E_b = \Gamma^a{}_{bc} E_a)$ に対応する（最後の $(c)$ が “微分方向”）。

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
[E_b, E_c] = - C^a{}_{bc} E_a.
$$

> 注意：多くの教科書では $([E_b,E_c]=+f^a{}_{bc}E_a)$ を採用する。
> 本プロジェクトでは **その $(f^a{}_{bc})$ に対して $(C^a{}_{bc}=-f^a{}_{bc})$** の規約を採用している。

### 3.3 runner 実装ルール（推奨）

* **C は手打ちしない**。可能なら runner 側で $de^a$ を明示し、係数比較で $C^a_{bc}$ を抽出して `self.data['structure_constants']=C` に入れる。
* 最低限、$C^a_{bc}$ が **b,c で反対称**になっていることを自動チェックすること。

## 4. 接続（スピン接続）とメトリック整合

接続 1-form を
$$
\omega^a{}_b = \Gamma^a{}_{bc} e^c
$$
で定義する。

メトリック整合（ローレンツ接続／直交接続）を仕様として固定：
$$
\omega_{ab} = -\omega_{ba}
\quad(\Leftrightarrow\quad
\Gamma_{abc} = -\Gamma_{bac})
$$
ただし $(\Gamma_{abc} = g_{ad}\Gamma^d{}_{bc})$。

## 5. Levi-Civita 接続（v3 engine の一般 Koszul 実装）

フレームが正規直交で、上記の構造定数規約を採用したとき、Levi-Civita 接続は v3 engine では次の**一般 Koszul 公式**で計算する：

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
T^a = de^a + \omega^a{}_b \wedge e^b,
\qquad
T^a = \frac12 T^a{}_{bc} e^b\wedge e^c.
$$

曲率 2-form：
$$
R^a{}_b = d\omega^a{}_b + \omega^a{}_c\wedge \omega^c{}_b,
\qquad
R^a{}_b = \frac12 R^a{}_{bcd} e^c\wedge e^d.
$$

## 7. 曲率成分の計算式（engine が実際に使う形）

v3 engine の $R^a_{bcd}$ は（現状）次の形を用いる：

$$
R^a{}_{bcd}
=
\Gamma^a{}_{ec}\Gamma^e{}_{bd}
-\Gamma^a{}_{ed}\Gamma^e{}_{bc}
-\Gamma^a{}_{be} C^e{}_{cd}.
$$

> 重要：一般にはここにフレーム方向微分項
> $(E_c(\Gamma^a{}_{bd}) - E_d(\Gamma^a{}_{bc}))$
> が現れるが、v3 engine ではそれを明示的に扱っていない。
> したがって runner は、左不変フレーム等により **$(\Gamma)$ がフレーム方向で定数扱い**になる設定を採用すること。

## 8. 必須セルフチェック（runner が満たすべき整合性）

runner は以下を満たすこと（落ちたら定義が engine と不整合）：

1. 構造定数の反対称：
   $$
   C^a{}_{bc} + C^a{}_{cb} = 0.
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


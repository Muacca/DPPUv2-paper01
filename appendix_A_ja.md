# 付録A：接続の成分式とコード入力の対応

## A.1 Levi–Civita接続

添字と構造係数は[第2章](section_02_ja.md)に従う。$S^3$の非零成分は式(A-01)である。

$$
\Gamma^{\mathrm{LC}}_{012}=-\Gamma^{\mathrm{LC}}_{021}
=\Gamma^{\mathrm{LC}}_{120}=-\Gamma^{\mathrm{LC}}_{102}
=\Gamma^{\mathrm{LC}}_{201}=-\Gamma^{\mathrm{LC}}_{210}
=\frac2r.
$$

$T^3$では全成分が零である。Nilの非零成分は式(A-02)となる。

$$
\Gamma^{\mathrm{LC}}_{012}=\Gamma^{\mathrm{LC}}_{021}=\frac1{2R},
\quad
\Gamma^{\mathrm{LC}}_{102}=\Gamma^{\mathrm{LC}}_{120}=-\frac1{2R},
\quad
\Gamma^{\mathrm{LC}}_{201}=-\Gamma^{\mathrm{LC}}_{210}=-\frac1{2R}.
$$

いずれも$C+\Gamma^{\mathrm{LC}}_{\cdot cb}-\Gamma^{\mathrm{LC}}_{\cdot bc}=0$を満たす。ここで空間スケールはframe位置によらず一定である。

## A.2 幾何学的torsionのcontortion

式(M2-03)へaxialとvector成分を代入すると、式(A-03)となる。

$$
K^{\mathrm{ax}}_{ijk}[T]=-\frac{\eta}{s}\varepsilon_{ijk},
\qquad
K^{\mathrm{vec}}_{abc}[T]
=-\frac13(\delta_{ac}V_b-\delta_{bc}V_a).
$$

axial成分の循環置換は偶置換なので、三項の和差は$T_{ijk}$そのものになる。vector成分では$\delta_{ab}V_c$の項が相殺され、上式を得る。両成分とも最初の二添字に反対称である。

## A.3 コード配列τとの対応

本稿に付属する計算コードの入力配列を$\tau_{abc}$と書く。この配列から式(A-04)により接続を構成する。

$$
K^{\mathrm{code}}_{abc}[\tau]
=\frac12(\tau_{abc}+\tau_{bca}-\tau_{cab}),\qquad
\Gamma^{\mathrm{code}}=\Gamma^{\mathrm{LC}}+K^{\mathrm{code}}.
$$

この接続からCartan方程式で再構成する幾何学的torsionは$T=-\tau$である。従って同じ接続の対応は式(A-05)となる。

$$
\tau=-T,\qquad
\eta_\tau=-\eta,\qquad V_\tau=-V,\qquad
K^{\mathrm{code}}[-T]=K[T].
$$

対応するパラメータを同時に変換すると、二つの表記で接続・曲率・作用は一致する。$K$だけを反転したり、$V$を据え置いて$\eta$だけを反転したりする変換ではない。

コード入力でのaxial contortionは$+\eta_\tau\varepsilon_{ijk}/s$、vector contortionは$+(\delta_{ac}(V_\tau)_b-\delta_{bc}(V_\tau)_a)/3$となる。幾何学的torsionの軸成分とvector成分は両方とも逆符号になる一方、$\eta^2,V^2,V\eta$は不変である。

## A.4 スケールと入力領域の変換

コード入力の区間$\eta_\tau\in[p,q]$は$\eta\in[-q,-p]$へ、$V_\tau\in[v_-,v_+]$は$V\in[-v_+,-v_-]$へ移る。図の軸を昇順にする場合には、ラベルだけでなく対応するデータ列の順序も反転する。

$S^3$の実半径と係数の辞書は式(A-06)である。

$$
a=\frac r2,\qquad
\frac{2\eta}{r}=\frac{2\chi}{a},\qquad
\chi=\frac{\eta}{2}=-\frac{\eta_\tau}{2}.
$$

$T^3$の周長表示では式(A-07)となる。

$$
\ell_i=2\pi R_i,\qquad
\frac{2\eta}{R_1}=\frac{2\zeta}{\ell_1},\qquad
\zeta=2\pi\eta=-2\pi\eta_\tau.
$$

Nilでは採用格子を固定し、共通スケール$R$を用いる。円方向の$L$は三幾何とも実円周長で、座標$u$の周期は1である。

## A.5 作用の比較に必要な区別

同時符号変換は同じ接続の表記変更である。NYの二形式係数、積分体積、作用variantの選択はそれぞれ別の定義であり、符号辞書で吸収しない。対応付ける両表記には同じ$n_{\mathrm{REE}}=\varepsilon R/2$と採用体積を用いる。

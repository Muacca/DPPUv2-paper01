# 7. 物理的解釈の範囲

## 7.1 縮約極小と全停留性

半径方向の厳密な極小であっても、省略した変分方向の勾配が零であるとは限らない。 $S^3$ の完全作用で $V=1,\eta=\theta=0,\kappa=L=1$ 、 $r=2\sqrt3$ を取ると、採用体積のもとで式(M7-01)となる。

$$
\begin{aligned}
\partial_r U&=0,&\partial_r^2U&=\sqrt3\thinspace \pi^2>0,\\
\partial_V U&=4\sqrt3\thinspace \pi^2\ne0,
&\partial_L U&=-4\sqrt3\thinspace \pi^2\ne0.
\end{aligned}
$$

$V$ と $L$ を固定した半径問題としては極小だが、それらの方向を許した停留点ではない。この例は、半径の停留条件だけではtorsion係数と円周長に関する停留条件を満たさないことを示す。

$T^3$ の異方変分にも同じ注意が必要である。固定形状経路の二階微分と、 $`(R_1,R_2,R_3)`$ を独立に動かしたHessianは異なる。[付録E](appendix_E_ja.md)では $`H_{22}=H_{33}=0`$ を明示し、半径方向の結果から全Hessian正定を結論しない。

## 7.2 S³×S¹の全真空に関する制約

ここでは、標準EC作用、定数係数の完全NY項、物質および外部spin sourceなし、宇宙定数なし、境界のない閉じた実Euclidean $S^3\times S^1$ 、滑らかで非退化なcoframe、独立な計量適合接続を仮定する。

完全NYは全微分であり、閉多様体上で定数係数を掛けてもbulk方程式を変えない。非退化なcoframeに対する標準ECの接続方程式はtorsion零を要求し、coframe方程式はRicci-flat条件になる。これらのbulk変分とNY恒等式については、[Corichi・Rubalcava-García・Vukašinac, *Actions, topological terms and boundaries in first order gravity: A review*, §IV.A–B、式(4.10)–(4.12)、(4.20)](https://arxiv.org/abs/1604.07764)を参照する。以下では、上記の実Euclidean・閉多様体の仮定のもとで幾何学的帰結を導く。

$S^3\times S^1$ のEuler数は零である。四次元の閉じたRiemannian多様体に対するGauss–Bonnet公式では、Ricci-flat条件のもとでEuler数がRiemann曲率ノルム二乗の積分の正定数倍となる。従って、その非負な被積分関数は至る所で零となり、計量は平坦でなければならない。

Bieberbachの定理により、コンパクト平坦四多様体は平坦な $T^4$ による有限被覆を持つ。その基本群には有限指数の $\mathbb Z^4$ 部分群が必要である（[LeBrun, *Einstein metrics, complex surfaces, and symplectic 4-manifolds*, Lemma 1の証明](https://doi.org/10.1017/S0305004109002345)も参照）。ところが $`\pi_1(S^3\times S^1)=\mathbb Z`$ の部分群は巡回群であり、 $\mathbb Z^4$ にはなれない。この矛盾から、上記の仮定を全て満たす真空解は存在しない。これは有限個のansatz点を調べた経験則ではなく、同じ作用・大域条件に対する制約である。

この結論は、 $S^3\times S^1$ をoff-shell背景として用いることを禁止しない。また、宇宙定数、source、境界、複素計量、非定数NY結合、高次曲率項を含む別の問題の不存在を証明するものでもない。それらの採用や解の探索は本稿の範囲に含めない。

## 7.3 NYと他の不変量

対象ansatzでの $`n_{\mathrm{TT}}=n_{\mathrm{REE}}`$ は、完全NY密度が零になる局所恒等式である。曲率側の $`\varepsilon^{abcd}R_{abcd}`$ は曲率に一次であり、曲率に二次のPontryagin形式 $`R^a{}_b\wedge R^b{}_a`$ とは異なる。torsion零の切片で曲率側の縮約が零になることを、全てのtorsion接続に対する消失へ拡張してはならない。

同様に、閉空間で完全NY積分が零であることから、個別のTT・REE密度が零であるとはいえない。ここで診断用作用が $\theta$ 依存性を持つことは、完全作用にNY由来のbulk効果が存在する根拠にはならない。

## 7.4 質量・崩壊・幾何間比較

$U=-S$ の符号やその極小から、経路積分を支配するsaddleを直ちに同定しない。積分路・重み・許容する変分を指定する必要がある。また、 $\Delta U$ を物理的障壁エネルギーや崩壊率へ変換するには、動的変数の運動項、制約、Euclidean境界条件および必要ならbounce解が必要である。type-Iの作用差だけではWKB指数も定まらない。

大域体積を固定したことで積分量の幾何学的正規化は定まったが、この物理的辞書は未確定である。本稿で比較するのは、同じ定義に基づく縮約関数とその極値構造である。最小の作用値を持つ幾何が物理的に選ばれるとするには、異なるトポロジーを比較する枠組みと共通の拘束条件が別途必要となる。

$S^2\times\mathbb R$ 等を対象へ広げる場合も、非コンパクト体積、境界、coframeと変分の指定が必要である。「 $S^2$ の接束が自明でない」ことだけから積空間のcoframeの不存在を結論することはできない。本稿はその拡張計算を行わない。

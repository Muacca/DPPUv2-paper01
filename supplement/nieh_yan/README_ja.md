# NY恒等式とT/τ対応の独立計算ノート

[English](README.md) · [計算ノートを開く](nieh_yan_identity.ipynb)

Paper 01の局所Nieh–Yan式と、幾何学的torsion $T$ ・コード入力 $\tau$ の対応を導く任意の補足資料です。定義・主張・論証は論文の第2章・付録A–Bで完結しており、このノートではその代数計算を実行しながら追えます。

## 計算内容

ノートには日英の説明を併記し、厳密な記号計算を用います。 $S^3$ ・ $T^3$ ・ $\mathrm{Nil}^3$ のcoframeを構成し、Cartan方程式から計量適合接続を解いたうえで、torsionと曲率を外微分形式として計算します。その後、次式の両辺を評価します。

$$
d(e^a\wedge T_a)=T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab}.
$$

成分展開からTTの係数 $1/4$ とREEの係数 $1/2$ が得られます。静的一様ansatzでは、三つの幾何とも

$$
n_{\mathrm{TT}}=n_{\mathrm{REE}}=-\frac{4V\eta}{s},
\qquad n_{\mathrm{FULL}}=0
$$

となります。 $s$ は球面では $r$ 、トーラスのansatzでは $`R_1`$ 、Nilでは $R$ で、球面の実半径は $r/2$ です。REEは正符号の曲率項、FULLはTTからREEを引いたものを表します。直接求めた接続を二つのcontortion式とも照合し、

$$
\tau=-T,\qquad V_\tau=-V,\qquad \eta_\tau=-\eta
$$

の対応を示します。入力領域もパラメータとともに反転させる必要があるため、区間端点の対応を説明しています。さらに混合torsionの厳密な数値例と、axialのみ・vectorのみの極限を示します。

## 実行方法

Python 3、SymPy、Jupyter Notebookを使用します。リポジトリ直下で実行してください。

```powershell
python -m pip install sympy notebook
jupyter notebook supplement/nieh_yan/nieh_yan_identity.ipynb
```

**Restart Kernel and Run All Cells**を選択します。配布ノートにはセルの実行結果を保存していません。密度行列の行はS3・T3・Nil3、列はTT・REE・FULL・ $`d(e^a\wedge T_a)/\mathrm{vol}_4`$ です。残差表示の空辞書は零形式を表し、成分表示との比較行列は零になります。

計算に使うのはSymPyとPythonの標準ライブラリだけです。データファイルの読込みや論文の計算エンジンの呼出しはなく、ノートを別のディレクトリへコピーしても実行できます。

## 論文との対応

- [英語論文](../../DPPUv2-paper01.pdf)：第2章・付録A–Bに作用、接続規約、成分係数、torsion辞書を記載しています。
- [日本語版第2章](../../section_02_ja.md)・[付録A](../../appendix_A_ja.md)・[付録B](../../appendix_B_ja.md)：対応する定義・導出です。
- [第3章](../../section_03_ja.md)・[規約文書](../../script/docs/CONVENTIONS_ja.md)：coframeの尺度、積分領域、大域体積を記載しています。

パラメータはframe位置によらない定数として扱います。局所微分形式の計算であり、FULLが点ごとに零になることはこのansatzの性質です。一般のNY恒等式と区別してください。大域正規化は論文で別途定めており、このノートでは相走査や全場方程式の求解を行いません。

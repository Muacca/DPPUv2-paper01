# 付録F：コード・データと再現手順

## F.1 コードと資料の提供先

コードと資料は[Muacca/DPPUv2-paper01](https://github.com/Muacca/DPPUv2-paper01)で提供する。リポジトリに付されたDOIは[10.5281/zenodo.18213677](https://doi.org/10.5281/zenodo.18213677)である。資料を引用する際には、使用した版とファイルを示す。

記号計算コアは[engine](script/DPPUv2_engine_core_v3.py)、幾何ごとの入口は[S³ runner](script/DPPUv2_runner_S3S1_v3.py)、[T³ runner](script/DPPUv2_runner_T3S1_v3.py)、[Nil³ runner](script/DPPUv2_runner_Nil3S1_v3.py)である。[パラメータ走査スクリプト](script/DPPUv2_parameter_scan_v3.py)は縮約ポテンシャルを評価する。仕様は[scriptの案内](script/README_ja.md)と[成分・大域規約](script/docs/CONVENTIONS_ja.md)、依存関係は[script/requirements.txt](script/requirements.txt)に示す。

## F.2 環境と計算例

以下はリポジトリ直下をカレントディレクトリとした実行例である。仮想環境を作る場合は`python -m venv .venv`を実行し、利用するshellで有効化する。依存パッケージは`python -m pip install -r script/requirements.txt`で導入する。

記号計算の例：

```powershell
python -B script/DPPUv2_runner_S3S1_v3.py --mode MX --ny-variant FULL --log-file new_s3_full.log
```

小さな走査の例：

```powershell
python -B script/DPPUv2_parameter_scan_v3.py --topologies S3 --ny-variants FULL --V-min 1 --V-max 4 --V-points 2 --eta-min -3 --eta-max 3 --eta-points 7 --theta-min 0 --theta-max 1 --theta-points 2 --workers 1 --output-dir new_scan
```

これは $`V_\tau>0`$ の入力例であり、本文では $V<0$ に対応する。ログと出力の保存先は実行ごとに分け、入力と結果を対応付ける。

## F.3 計算条件の記録

コードの版、作用variant、 $T/\tau$ の入力規約、体積規約、変分経路、入力の点数と端点、分類器の境界処理を記録する。 $S^3$ の $r$ と実半径 $a$ 、 $T^3$ の $`R_i`$ と周長 $`\ell_i`$ を分ける。Nilでは格子の面同定も含めて空間を特定する。

CSVの8列と入力域は[付録C](appendix_C_ja.md)、表示量の定義は[付録D](appendix_D_ja.md)に示した。notebookを用いる場合、出力には各セルを実行した時点の値が格納される。計算の再現には、使用する環境を指定し、kernelを再起動したうえで必要な計算セルを順に実行する。

## F.4 結果の記録と引用

計算結果には、幾何、作用variant、パラメータ値、コードの版を添える。CSVの数値、ポテンシャル曲線、相分類が同じ入力を指すことを確認し、式による解析結果と入力格子上の数値結果を区別して記載する。

ライセンスは[LICENSE](LICENSE)に従う。計算手順・入力条件・出力を一組として提示することで、同じ条件での再現と比較が可能になる。

# DPPUv2 Paper 01：アインシュタイン–カルタン重力における縮約ポテンシャルの幾何依存分類

[English README](README.md)

本リポジトリは論文の補足資料として、本文、原稿ソース、計算コード、掲載図表のデータを提供します。本論文はEinstein–Cartan重力において、$S^3\times S^1$、$T^3\times S^1$、コンパクト商空間$\mathrm{Nil}^3\times S^1$上の静的・一様な幾何に対する縮約ポテンシャルを比較します。

対象ansatzでは完全Nieh–Yan（NY）密度が零となり、FULLの縮約ポテンシャルはNY結合に依存しません。TTと正符号REEは別の診断用作用として扱います。相分類は指定した変分経路上のoff-shell極値の分類であり、半径方向の極小だけで全真空解、物理的な質量や崩壊率が確立されるわけではありません。

## 論文の版

**v3 — 改訂日：2026-10-25。** 本文・図表・付属データはv3に対応しています。v2からの作用、torsionの表記、大域正規化、結論の変更は[改訂履歴](CHANGELOG.md)にまとめています。

リポジトリのDOI：[10.5281/zenodo.18213677](https://doi.org/10.5281/zenodo.18213677)。引用や結果の比較では、使用した原稿の版とファイルを明記してください。

[訂正追跡（Issue #4）](https://github.com/Muacca/DPPUv2-paper01/issues/4)

## 論文・データを読む

| 資料 | 内容 |
|---|---|
| [英語PDF](DPPUv2-paper01.pdf) | 本文・図表・付録を含むv3論文 |
| [日本語Markdown](section_01_ja.md) | v3論文。序論から各章・付録へ移動できます |
| [掲載図表のデータ](data/paper_v3/README.md) | 数値データ、パラメータ、入力域、図一覧 |
| [LaTeXソース](LaTeX/main.tex) | メインファイル、各章、付録 |
| [計算コード](script/) | 記号計算エンジン、各幾何のrunner、パラメータ走査、可視化ツール |
| [規約](script/docs/CONVENTIONS_ja.md) | 接続・torsion・添字と大域幾何の規約 |
| [NY独立計算ノート](supplement/nieh_yan/README_ja.md) | NY恒等式・成分係数・T/τ対応を外積計算で追う任意の補足資料 |

![完全NY作用の相分類](LaTeX/figures/v3/fig05_full_phases.png)

この相図の軸は幾何学的torsionのパラメータです。定義と変分経路は本文を参照してください。

## 計算コードを使う

以下はリポジトリ直下をカレントディレクトリとしたコマンドです。[script/requirements.txt](script/requirements.txt)の依存パッケージを導入したPython環境を使用してください。

```powershell
python -m pip install -r script/requirements.txt
python -B script/DPPUv2_parameter_scan_v3.py --help
```

[付録F](appendix_F_ja.md)に示した小さな走査の例：

```powershell
python -B script/DPPUv2_parameter_scan_v3.py --topologies S3 --ny-variants FULL --V-min 1 --V-max 4 --V-points 2 --eta-min -3 --eta-max 3 --eta-points 7 --theta-min 0 --theta-max 1 --theta-points 2 --workers 1 --output-dir new_scan
```

本文と掲載データは幾何学的torsion $T$を使います。コード入力は$\tau=-T$であり、$V_\tau=-V$、$\eta_\tau=-\eta$です。したがって、この例の正のV入力は本文の負のVに対応します。比較時はtorsion規約、大域体積、変分経路を揃えてください。パラメータと入力域の対応は[付録A](appendix_A_ja.md)・[付録C](appendix_C_ja.md)、再現手順は付録Fに示しています。

CSVは指定ディレクトリに`dppu_scan_all_{timestamp}.csv`および`dppu_scan_{topology}_{ny_variant}_{timestamp}.csv`として保存されます。走査出力と掲載図表のデータでは形式が異なるため、付録Cと[データ案内](data/paper_v3/README.md)を参照してください。

計算ごとに新しい出力・checkpointディレクトリを使い、コードの版を変えた場合はnotebookやviewerのkernelを再起動してください。notebookには実行済み出力を含めていません。用途と実行条件は[スクリプト案内](script/README_ja.md)を参照してください。v3の図は[LaTeX/figures/v3](LaTeX/figures/v3/)、表は[LaTeX/tables/v3](LaTeX/tables/v3/)、対応するデータは[data/paper_v3](data/paper_v3/)にあります。改訂前の結果とその解釈上の制限は[改訂履歴](CHANGELOG.md#earlier-results)を参照してください。

## PDFをビルドする

`pdflatex`と[main.tex](LaTeX/main.tex)で指定するパッケージを備えたLaTeX環境を使用します。リポジトリ直下から次を実行してください。

```powershell
Push-Location LaTeX
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
Pop-Location
```

複数回のコンパイルで相互参照と目次を解決します。生成先は`LaTeX/main.pdf`です。エラーや未解決の参照は`LaTeX/main.log`で確認してください。このコマンドはリポジトリ直下の配布用`DPPUv2-paper01.pdf`を置き換えません。

## 著者・ライセンス

Muacca。利用条件は[LICENSE](LICENSE)を参照してください。

# DPPUv2 Phase 1 - Nieh-Yan項を含むEinstein-Cartan理論

このディレクトリには、DPPUv2（Differential Geometric Phase Portrait with Uniaxial Torsion, version 2）Phase 1の実装が含まれています。Nieh-Yan位相項を含むEinstein-Cartan理論の記号計算と数値解析を行います。

⇒ [English version](README.md)

## 概要

DPPUv2 Phase 1は以下の機能を提供します：
- ねじれを含むEinstein-Cartan場の方程式の記号計算
- 3種類の時空トポロジーの解析：S³×S¹、T³×S¹、Nil³×S¹
- 3種類のねじれアンザッツモード：軸性（AX）、ベクトル・トレース（VT）、混合（MX）
- 完全Nieh-Yan項（FULL）と診断用作用（TT・正符号REE）の比較
- 相図の生成と可視化
- 全パラメータ組み合わせの一括実行

## コアスクリプト

### エンジンとトポロジーランナー

- **DPPUv2_engine_core_v3.py** - Einstein-Cartan理論を実装するコア計算エンジン
  - ねじれアンザッツ計算の基礎フレームワークを提供
  - Riemannテンソルの反対称性の厳密な検証を実装
  - 3つのモード（AX/VT/MX）と3つの作用選択（FULLと診断用のTT/REE）をサポート
  - 型安全なenumベースのモード管理

- **DPPUv2_runner_S3S1_v3.py** - S³×S¹トポロジーランナー
  - 3次元球面×円トポロジー
  - 構造定数：C^i_jk = (4/r)ε_ijk（SU(2)リー代数）
  - 曲率：R_LC = 24/r²

- **DPPUv2_runner_T3S1_v3.py** - T³×S¹トポロジーランナー
  - 3次元トーラス×円トポロジー
  - 平坦な構造定数：C^i_jk = 0
  - 曲率：R_LC = 0

- **DPPUv2_runner_Nil3S1_v3.py** - Nil³×S¹トポロジーランナー
  - Nil多様体×円トポロジー
  - 非零の構造係数： $`C^2{}_{01}=-1/R`$ 、 $`C^2{}_{10}=1/R`$ （空間添字 $0,1,2$ 、スケール $R>0$ ）
  - Levi–Civitaスカラー曲率： $`\mathcal R_{\mathrm{LC}}=-\frac{1}{2R^2}`$ （torsionを含まない接続）

### 解析ツール

- **DPPUv2_parameter_scan_v3.py** - 相図データ生成
  - (V, η, θ_NY)パラメータ空間を系統的にスキャン
  - 安定性分類（type-I/II/III）を含むCSVファイルを出力
  - 並列計算のためのマルチプロセッシングをサポート

- **DPPUv2_visualize_phasemap_v3.py** - 相図可視化
  - CSVデータから相図画像を生成
  - カラーコード化された安定性マップを作成
  - θ_NY変化に対するマルチフレーム可視化をサポート

- **DPPUv2_visualize_phasematrix_v3.py** - 相行列可視化
  - 複数パラメータにわたる包括的な行列プロットを作成
  - トポロジーとバリアント間の比較解析

## ディレクトリ構造

### `bin/` - 一括実行スクリプト

全ケースのrunner を一括実行するためのスクリプト：
- `run_all.bat` - 自動テスト用Windowsバッチファイル
- `run_all.sh` - 自動テスト用Linux/macOS bashスクリプト
- `README.md` - 一括実行のドキュメント（英語）
- `README_ja.md` - 一括実行のドキュメント（日本語）

詳細な使用方法は [bin/README_ja.md](bin/README_ja.md) を参照してください。

### `docs/` - ドキュメント

技術ドキュメントと規約：
- [DPPUv2 Engine CONVENTIONS](docs/CONVENTIONS_ja.md) - DPPUv2 エンジンの規約と仕様
- [DPPUv2 SymPy guideline](docs/SymPy_guideline_ja.md) - SymPy使用ガイドラインとベストプラクティス

### 掲載データとnotebook

[v3データ案内](../data/paper_v3/README.md)に、掲載図表・パラメータ・数値データをまとめています。このデータは内容別の列と幾何学的torsionのパラメータを使っており、CSV可視化スクリプトが読み込むscannerの8列CSVとは形式が異なります。

- [対話型可視化](DPPUv2_visualize_notebook_v3.ipynb)：付属viewerを起動します。計算用の依存ライブラリに加え、Jupyter・ipywidgets・ipymplが必要です。
- [ポテンシャルの探索](DPPUv2_Paper_Figures_v3.ipynb)：付属エンジンを使って項別分解と完全NY作用・診断用作用の比較を行います。記載した膨張経路ではFULLのtype-Iフィルタは空集合を返します。論文掲載の9図を生成するnotebookではありません。

いずれも実行済み出力を含めていません。この`script/`ディレクトリをカレントディレクトリとしてJupyterを起動し、新しいkernelで実行してください。図のnotebookは描画セルの実行時にカレントディレクトリへPNG/PDFを出力します。notebookの変数`V`・`eta`はコード入力 $`V_\tau,\eta_\tau`$ で、本文では $`V=-V_\tau`$ 、 $`\eta=-\eta_\tau`$ です。半径・体積の規約は[CONVENTIONS](docs/CONVENTIONS_ja.md)、改訂前の出力は[改訂履歴](../CHANGELOG.md#earlier-results)を参照してください。

## クイックスタート

以下のコマンドは、この`script/`ディレクトリをカレントディレクトリとして実行します。

### 単一トポロジーの実行

選択したモードとNYバリアントで特定のトポロジーを実行：

```bash
# S³×S¹トポロジー、混合モード、FULLバリアント
python DPPUv2_runner_S3S1_v3.py --mode MX --ny-variant FULL

# T³×S¹トポロジー、軸性モード、TTバリアント
python DPPUv2_runner_T3S1_v3.py --mode AX --ny-variant TT

# Nil³×S¹トポロジー、ベクトル・トレースモード、REEバリアント
python DPPUv2_runner_Nil3S1_v3.py --mode VT --ny-variant REE
```

### 一括実行（全組み合わせ）

**Windows:**
```cmd
bin\run_all.bat          # フルモード：27ケース
bin\run_all.bat quick    # クイックモード：3ケース
```

**Linux/macOS:**
```bash
chmod +x bin/run_all.sh
bin/run_all.sh         # フルモード：27ケース
bin/run_all.sh quick   # クイックモード：3ケース
```

### パラメータスキャンと相図生成

```bash
# 相図データを生成
python -B DPPUv2_parameter_scan_v3.py --topologies S3 --ny-variants FULL --output-dir output

# YYYYMMDD_HHMMSSを生成されたCSVの日時に置き換えて可視化
python -B DPPUv2_visualize_phasemap_v3.py output/dppu_scan_S3_FULL_YYYYMMDD_HHMMSS.csv
```

## コマンドラインオプション

### トポロジーランナー

```
--mode {AX,VT,MX}           ねじれアンザッツモード
                            AX: 軸性のみ（S^μ ≠ 0, T_μ = 0）
                            VT: ベクトル・トレースのみ（T_μ ≠ 0, S^μ = 0）
                            MX: 混合（両方とも非ゼロ）

--ny-variant {TT,REE,FULL}  作用の選択（定義は下記）

--log-file PATH             出力ログファイルのパス
--checkpoint-dir PATH       計算チェックポイント用ディレクトリ
```

EC項に加える項を次のように選びます。各項は四次元で積分します。

- `TT`：torsion二形式の積 $`T^a\wedge T_a`$ に $+\theta$ を掛けた診断用作用。
- `REE`：曲率側 $`e^a\wedge e^b\wedge R_{ab}`$ に $+\theta$ を掛けた診断用作用。
- `FULL`：両者の差 $`T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab}`$ に $\theta$ を掛けた完全Nieh-Yan項。

ここで $T$ はCartan方程式で定義する幾何学的torsionです。TTはvector trace成分だけを指す名称ではありません。接続・入力配列との対応は[規約文書](docs/CONVENTIONS_ja.md)、作用の定義は[論文第2章](../section_02_ja.md)を参照してください。

### パラメータスキャン

```
--topologies {S3,T3,Nil3} [...]    対象幾何（複数指定可）
--ny-variants {FULL,TT,REE} [...]  作用の選択（複数指定可）
--output-dir DIR                  出力CSVを保存するディレクトリ
--workers N                       並列ワーカー数
```

入力範囲・点数のオプションは`python -B DPPUv2_parameter_scan_v3.py --help`で確認できます。出力は指定ディレクトリの`dppu_scan_all_{timestamp}.csv`および`dppu_scan_{topology}_{ny_variant}_{timestamp}.csv`です。入力の意味・既定値と少数点の実行例は[付録C](../appendix_C_ja.md)・[付録F](../appendix_F_ja.md)に示します。

## 出力ファイル

### ランナー出力
- **LOGファイル**: 記号式を含む詳細な計算ログ
- **SUCCESSマーカー**: 計算完了の成功を示す
- **チェックポイントファイル**: 中間計算状態（有効化時）

### 一括実行出力
- **summary_*.log**: 成功/失敗統計を含む全体実行サマリー
- **results_*.csv**: CSV形式の結果テーブル
- **個別ケースログ**: 各テストケースの詳細ログ

### 相図出力
- **CSVデータファイル**: パラメータ空間スキャン結果
- **PNG画像**: 相図可視化

## 依存関係

- Python 3.7+
- SymPy: 記号数学
- NumPy: 数値計算
- SciPy: 最適化アルゴリズム
- Matplotlib: 可視化（相図用）
- Pandas: データ処理（相図用）

## 参考文献

- Chandia & Zanelli (1997): Phys. Rev. D 55, 7580 [hep-th/9702025]
- Hehl et al. (1976): Rev. Mod. Phys. 48, 393
- Nieh & Yan (1982): J. Math. Phys. 23, 373

## バージョン履歴

- **v3.0** (2025-12-14): 論文のための大規模リファクタリング
  - モード名：AX/VT/MX
  - NYバリアント：TT/REE/FULL
  - Riemann反対称性の厳密検証（3段階チェック）
  - enumベースの型安全なモード管理
  - インフラストラクチャの分離（logger/checkpoint はオプション）

## 作者

Muacca

## ライセンス

リポジトリルートのLICENSEファイルを参照してください。

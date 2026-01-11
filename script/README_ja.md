# DPPUv2 Phase 1 - Nieh-Yan項を含むEinstein-Cartan理論

このディレクトリには、DPPUv2（Differential Geometric Phase Portrait with Uniaxial Torsion, version 2）Phase 1の実装が含まれています。Nieh-Yan位相項を含むEinstein-Cartan理論の記号計算と数値解析を行います。

⇒ [English version](README.md)

## 概要

DPPUv2 Phase 1は以下の機能を提供します：
- ねじれを含むEinstein-Cartan場の方程式の記号計算
- 3種類の時空トポロジーの解析：S³×S¹、T³×S¹、Nil³×S¹
- 3種類のねじれアンザッツモード：軸性（AX）、ベクトル・トレース（VT）、混合（MX）
- 3種類のNieh-Yan結合バリアント：TT、REE、FULL
- 相図の生成と可視化
- 全パラメータ組み合わせの一括実行

## コアスクリプト

### エンジンとトポロジーランナー

- **DPPUv2_engine_core_v3.py** - Einstein-Cartan理論を実装するコア計算エンジン
  - ねじれアンザッツ計算の基礎フレームワークを提供
  - Riemannテンソルの反対称性の厳密な検証を実装
  - 3つのモード（AX/VT/MX）と3つのNYバリアント（TT/REE/FULL）をサポート
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
  - 構造定数：C¹_23 = -C¹_32 = 4/r（Heisenberg代数）
  - 曲率：R_LC = 0

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
- [DPPUv2_engine_core_v3_CONVENTIONS](docs/DPPUv2_engine_core_v3_CONVENTIONS_ja.md) - エンジンコアの規約と仕様
- [DPPUv2_SymPy_guideline_v1](docs/DPPUv2_SymPy_guideline_v1_ja.md) - SymPy使用ガイドラインとベストプラクティス

### `sample_result/` - run_all/scan 実行結果サンプル

- `sample_result-run_all_20251221.zip`: 全ケースを一括実行した結果のサンプル
    - `DPPUv2_parameter_scan_v3` はこの結果をもとに作成
- `sample_result-scan_20251221.zip`: `DPPUv2_parameter_scan_v3` の実行結果サンプル

## クイックスタート

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
python DPPUv2_parameter_scan_v3.py --topology S3S1 --ny-variant FULL

# 結果を可視化
python DPPUv2_visualize_phasemap_v3.py results_data.csv
```

## コマンドラインオプション

### トポロジーランナー

```
--mode {AX,VT,MX}           ねじれアンザッツモード
                            AX: 軸性のみ（S^μ ≠ 0, T_μ = 0）
                            VT: ベクトル・トレースのみ（T_μ ≠ 0, S^μ = 0）
                            MX: 混合（両方とも非ゼロ）

--ny-variant {TT,REE,FULL}  Nieh-Yan結合バリアント
                            TT: ねじれ・トレースのみ
                            REE: Riemannのみ
                            FULL: 両項

--log-file PATH             出力ログファイルのパス
--checkpoint-dir PATH       計算チェックポイント用ディレクトリ
```

### パラメータスキャン

```
--topology {S3S1,T3S1,Nil3S1}  時空トポロジー
--ny-variant {TT,REE,FULL}     Nieh-Yanバリアント
--output PATH                  出力CSVファイルパス
--workers N                    並列ワーカー数
```

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

## 著者

Muacca

## ライセンス

リポジトリ直下のLICENSEファイルを参照してください。

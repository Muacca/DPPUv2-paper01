# DPPUv2 v3.0 一括実行スクリプト

全モード・全トポロジーを一括実行し、結果を出力するスクリプトです。

⇒ [English version](README.md)

## ファイル

- **run_all.bat** - Windows バッチファイル版
- **run_all.sh** - Linux/macOS bash スクリプト版

## 使用方法

### Windows バッチファイル版

```cmd
# 全テスト実行（27ケース）
bin\run_all.bat

# クイックモード（3ケース）
bin\run_all.bat quick
```

### Linux/macOS bash スクリプト版

```bash
# 実行権限の付与（初回のみ）
chmod +x bin/run_all.sh

# 全テスト実行（27ケース）
bin/run_all.sh

# クイックモード（3ケース）
bin/run_all.sh quick
```

## テストケース

### フルモード（27ケース）
- トポロジー: S3xS1, T3xS1, Nil3xS1（3種類）
- モード: AX, VT, MX（3種類）
- バリアント: TT, REE, FULL（3種類）
- 合計: 3 × 3 × 3 = **27ケース**

### クイックモード（3ケース）
- トポロジー: S3xS1, T3xS1, Nil3xS1（3種類）
- モード: MX（1種類）
- バリアント:  FULL（1種類）
- 合計: 3 × 1 × 1 = **3ケース**

## 出力ファイル

実行すると `run_results_YYYYMMDD_HHMMSS/` ディレクトリが作成され、以下のファイルが生成されます:

```
run_results_20251214_131718/
├── summary_20251214_131718.log       # サマリーログ
├── results_20251214_131718.csv       # CSV形式の結果
├── S3S1_AX_TT_20251214_131718.log    # 各ケースのログ
├── S3S1_AX_REE_20251214_131718.log
├── S3S1_AX_FULL_20251214_131718.log
└── ...
```

### summary_*.log の内容

- 実行時刻
- 全ケース数、成功/失敗数、成功率
- 各ケースの実行結果（PASS/FAIL）
- 失敗ケースの詳細（ある場合）
- 実行時間統計（トポロジー別の平均/最小/最大）

### results_*.csv の内容

CSV形式で以下の情報を記録:

| Topology | Mode | Variant | Status | Runtime(s) | LogFile | ErrorMessage |
|----------|------|---------|--------|------------|---------|--------------|
| S3S1     | AX   | TT      | PASS   | 0.58       | S3S1_AX_TT_*.log | "" |
| ...      | ...  | ...     | ...    | ...        | ...     | ... |

## 終了コード

- **0**: 全テスト成功
- **1**: 1件以上のテスト失敗


## トラブルシューティング

## 関連ファイル

- `DPPUv2_runner_S3S1_v3.py`   - S3×S1 runner
- `DPPUv2_runner_T3S1_v3.py`   - T3×S1 runner  
- `DPPUv2_runner_Nil3S1_v3.py` - Nil3×S1 runner
- `DPPUv2_engine_core_v3.py`   - 共通エンジン



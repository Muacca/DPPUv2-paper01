# DPPUv2 v3.0 Batch Execution Script

A script that executes all modes and all topologies in batch and outputs the results.

⇒ [日本語版](README_ja.md)

## Files

- **run_all.bat** - Windows batch file version
- **run_all.sh** - Linux/macOS bash script version

## Usage

### Windows Batch File Version

```cmd
# Run all tests (27 cases)
bin\run_all.bat

# Quick mode (3 cases)
bin\run_all.bat quick
```

### Linux/macOS Bash Script Version

```bash
# Grant execution permission (first time only)
chmod +x bin/run_all.sh

# Run all tests (27 cases)
bin/run_all.sh

# Quick mode (3 cases)
bin/run_all.sh quick
```

## Test Cases

### Full Mode (27 cases)
- Topology: S3xS1, T3xS1, Nil3xS1 (3 types)
- Mode: AX, VT, MX (3 types)
- Variant: TT, REE, FULL (3 types)
- Total: 3 × 3 × 3 = **27 cases**

### Quick Mode (3 cases)
- Topology: S3xS1, T3xS1, Nil3xS1 (3 types)
- Mode: MX (1 type)
- Variant: FULL (1 type)
- Total: 3 × 1 × 1 = **3 cases**

## Output Files

When executed, a `run_results_YYYYMMDD_HHMMSS/` directory is created with the following files:

```
run_results_20251214_131718/
├── summary_20251214_131718.log       # Summary log
├── results_20251214_131718.csv       # Results in CSV format
├── S3S1_AX_TT_20251214_131718.log    # Log for each case
├── S3S1_AX_REE_20251214_131718.log
├── S3S1_AX_FULL_20251214_131718.log
└── ...
```

### summary_*.log Contents

- Execution timestamp
- Total number of cases, success/failure counts, success rate
- Execution results for each case (PASS/FAIL)
- Details of failed cases (if any)
- Runtime statistics (average/min/max by topology)

### results_*.csv Contents

Records the following information in CSV format:

| Topology | Mode | Variant | Status | Runtime(s) | LogFile | ErrorMessage |
|----------|------|---------|--------|------------|---------|--------------|
| S3S1     | AX   | TT      | PASS   | 0.58       | S3S1_AX_TT_*.log | "" |
| ...      | ...  | ...     | ...    | ...        | ...     | ... |

## Exit Codes

- **0**: All tests passed
- **1**: One or more tests failed

This script is suitable for use in CI/CD pipelines.

## Troubleshooting

## Related Files

- `DPPUv2_runner_S3S1_v3.py`   - S3×S1 runner
- `DPPUv2_runner_T3S1_v3.py`   - T3×S1 runner  
- `DPPUv2_runner_Nil3S1_v3.py` - Nil3×S1 runner
- `DPPUv2_engine_core_v3.py`   - Common engine

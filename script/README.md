# DPPUv2 Phase 1 - Einstein-Cartan Theory with Nieh-Yan Term

This directory contains the implementation of DPPUv2 (Differential Geometric Phase Portrait with Uniaxial Torsion, version 2) Phase 1, which performs symbolic computation and numerical analysis of Einstein-Cartan theory with the Nieh-Yan topological term.

⇒ [日本語版](README_ja.md)

## Overview

DPPUv2 Phase 1 provides tools for:
- Symbolic computation of Einstein-Cartan field equations with torsion
- Analysis of three different spacetime topologies: S³×S¹, T³×S¹, and Nil³×S¹
- Three torsion ansatz modes: Axial (AX), Vector-trace (VT), and Mixed (MX)
- Three Nieh-Yan coupling variants: TT, REE, and FULL
- Phase diagram generation and visualization
- Batch execution across all parameter combinations

## Core Scripts

### Engine and Topology Runners

- **DPPUv2_engine_core_v3.py** - Core computation engine implementing Einstein-Cartan theory
  - Provides base framework for torsion ansatz computations
  - Implements strict Riemann tensor antisymmetry verification
  - Supports three modes (AX/VT/MX) and three NY variants (TT/REE/FULL)
  - Type-safe enum-based mode management

- **DPPUv2_runner_S3S1_v3.py** - S³×S¹ topology runner
  - 3-sphere × circle topology
  - Structure constants: C^i_jk = (4/r)ε_ijk (SU(2) Lie algebra)
  - Curvature: R_LC = 24/r²

- **DPPUv2_runner_T3S1_v3.py** - T³×S¹ topology runner
  - 3-torus × circle topology
  - Flat structure constants: C^i_jk = 0
  - Curvature: R_LC = 0

- **DPPUv2_runner_Nil3S1_v3.py** - Nil³×S¹ topology runner
  - Nilmanifold × circle topology
  - Structure constants: C¹_23 = -C¹_32 = 4/r (Heisenberg algebra)
  - Curvature: R_LC = 0

### Analysis Tools

- **DPPUv2_parameter_scan_v3.py** - Phase diagram data generation
  - Systematically scans (V, η, θ_NY) parameter space
  - Outputs CSV files with stability classifications (type-I/II/III)
  - Supports multiprocessing for parallel computation
  
- **DPPUv2_visualize_phasemap_v3.py** - Phase diagram visualization
  - Generates phase diagram images from CSV data
  - Creates color-coded stability maps
  - Supports multi-frame visualization for θ_NY variations

- **DPPUv2_visualize_phasematrix_v3.py** - Phase matrix visualization
  - Creates comprehensive matrix plots across multiple parameters
  - Comparative analysis across topologies and variants

## Directory Structure

### `bin/` - Batch Execution Scripts

Contains scripts for running all cases runner in batch:
- `run_all.bat` - Windows batch file for automated testing
- `run_all.sh` - Linux/macOS bash script for automated testing
- `README.md` - Documentation for batch execution (English)
- `README_ja.md` - Documentation for batch execution (Japanese)

See [bin/README.md](bin/README.md) for detailed usage instructions.

### `docs/` - Documentation

Technical documentation and conventions:
- [DPPUv2_engine_core_v3_CONVENTIONS](docs/DPPUv2_engine_core_v3_CONVENTIONS.md) - Engine core conventions and specifications
- [DPPUv2_SymPy_guideline_v1](docs/DPPUv2_SymPy_guideline_v1.md) - SymPy usage guidelines and best practices

### `sample_result/` - Sample run_all/scan execution results

- `sample_result-run_all_20251221.zip`: Sample results of running all cases at once
    - `DPPUv2_parameter_scan_v3` was created based on these results
- `sample_result-scan_20251221.zip`:  Sample results of `DPPUv2_parameter_scan_v3`

## Quick Start

### Single Topology Execution

Run a specific topology with chosen mode and NY variant:

```bash
# S³×S¹ topology, Mixed mode, FULL variant
python DPPUv2_runner_S3S1_v3.py --mode MX --ny-variant FULL

# T³×S¹ topology, Axial mode, TT variant
python DPPUv2_runner_T3S1_v3.py --mode AX --ny-variant TT

# Nil³×S¹ topology, Vector-trace mode, REE variant
python DPPUv2_runner_Nil3S1_v3.py --mode VT --ny-variant REE
```

### Batch Execution (All Combinations)

**Windows:**
```cmd
bin\run_all.bat          # Full mode: 27 cases
bin\run_all.bat quick    # Quick mode: 3 cases
```

**Linux/macOS:**
```bash
chmod +x bin/run_all.sh
bin/run_all.sh         # Full mode: 27 cases
bin/run_all.sh quick   # Quick mode: 3 cases
```

### Parameter Scan and Phase Diagram

```bash
# Generate phase diagram data
python DPPUv2_parameter_scan_v3.py --topology S3S1 --ny-variant FULL

# Visualize results
python DPPUv2_visualize_phasemap_v3.py results_data.csv
```

## Command Line Options

### Topology Runners

```
--mode {AX,VT,MX}           Torsion ansatz mode
                            AX: Axial-only (S^μ ≠ 0, T_μ = 0)
                            VT: Vector-trace-only (T_μ ≠ 0, S^μ = 0)
                            MX: Mixed (both non-zero)

--ny-variant {TT,REE,FULL}  Nieh-Yan coupling variant
                            TT: Torsion-trace only
                            REE: Riemann only
                            FULL: Both terms

--log-file PATH             Output log file path
--checkpoint-dir PATH       Directory for computation checkpoints
```

### Parameter Scan

```
--topology {S3S1,T3S1,Nil3S1}  Spacetime topology
--ny-variant {TT,REE,FULL}     Nieh-Yan variant
--output PATH                  Output CSV file path
--workers N                    Number of parallel workers
```

## Output Files

### Runner Output
- **LOG files**: Detailed computation logs with symbolic expressions
- **SUCCESS markers**: Indicates successful computation completion
- **Checkpoint files**: Intermediate computation states (if enabled)

### Batch Execution Output
- **summary_*.log**: Overall execution summary with pass/fail statistics
- **results_*.csv**: CSV format results table
- **Individual case logs**: Detailed logs for each test case

### Phase Diagram Output
- **CSV data files**: Parameter space scan results
- **PNG images**: Phase diagram visualizations

## Dependencies

- Python 3.7+
- SymPy: Symbolic mathematics
- NumPy: Numerical computations
- SciPy: Optimization algorithms
- Matplotlib: Visualization (for phase diagrams)
- Pandas: Data handling (for phase diagrams)

## References

- Chandia & Zanelli (1997): Phys. Rev. D 55, 7580 [hep-th/9702025]
- Hehl et al. (1976): Rev. Mod. Phys. 48, 393
- Nieh & Yan (1982): J. Math. Phys. 23, 373

## Version History

- **v3.0** (2025-12-14): Major refactoring for publication
  - Mode names: AX/VT/MX
  - NY variants: TT/REE/FULL
  - Strict Riemann antisymmetry verification (3-level check)
  - Enum-based type-safe mode management
  - Infrastructure separation (logger/checkpoint optional)

## Author

Muacca

## License

See LICENSE file in the repository root.

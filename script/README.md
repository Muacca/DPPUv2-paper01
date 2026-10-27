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
  - Nonzero structure coefficients: $`C^2{}_{01}=-1/R`$, $`C^2{}_{10}=1/R`$ (spatial indices $0,1,2$; scale $R>0$)
  - Levi–Civita scalar curvature: $`\mathcal R_{\mathrm{LC}}=-\frac{1}{2R^2}`$ (torsion-free connection)

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
- [DPPUv2 Engine CONVENTIONS](docs/CONVENTIONS.md) - DPPUv2 Engine conventions and specifications
- [DPPUv2 SymPy guideline](docs/SymPy_guideline.md) - SymPy usage guidelines and best practices

### Paper Data and Notebooks

The [v3 data guide](../data/paper_v3/README.md) lists the published figures, tables, parameters and numerical datasets. These datasets use descriptive columns and geometric-torsion parameters; they are not the eight-column scanner CSV format expected by the CSV visualization scripts.

- [Interactive visualizer](DPPUv2_visualize_notebook_v3.ipynb): launches the accompanying viewer. Requires Jupyter, ipywidgets and ipympl in addition to the calculation dependencies.
- [Paper figures](DPPUv2_Paper_Figures_v3.ipynb): generates and displays Figures 1–9 of the v3 paper. By default it reads the supplied CSV datasets; an optional, separate step recomputes the same fixed grids from the reduced polynomials.

Both notebooks are supplied without evaluated outputs. Start Jupyter from this `script/` directory with a fresh kernel. The figure notebook also accepts the repository root as its working directory. It uses geometric-torsion parameters, while the interactive viewer and engine APIs use code-input parameters. The dictionary is $`V_\tau=-V,\eta_\tau=-\eta`$. See [CONVENTIONS](docs/CONVENTIONS.md) for the radius and volume conventions.

### Reproducing Figures 1–9

Use Python 3.12 and install `python -m pip install -r requirements-figures.txt` from this directory. This file pins NumPy 2.3.5, Matplotlib 3.10.7 and Pillow 12.0.0 for the figure environment; it does not replace the engine's `requirements.txt`. The renderer uses DejaVu Serif and saves PNG at 360 dpi plus SVG.

In the notebook, choose a new `OUTPUT_ROOT` and leave `RECOMPUTE = False` to draw from `data/paper_v3`. Set `RECOMPUTE = True` to first write fresh CSVs to `OUTPUT_ROOT/data`. Each numbered cell then saves and displays one figure in `OUTPUT_ROOT/figures`. The default run directory is `script/paper_v3_run`. Existing outputs are not overwritten; choose a new run directory or display an already generated PNG when revisiting a cell.

The same operations are available without Jupyter, from `script/`:

```powershell
# Draw directly from the supplied data.
python -B -m paper_figures render --output-dir new_figures

# Separately recompute the fixed v3 datasets, then draw from them.
python -B -m paper_figures recompute --output-dir new_data
python -B -m paper_figures render --data-dir new_data --output-dir new_figures_recomputed
```

[paper_figures/data.py](paper_figures/data.py) evaluates the fixed figure inputs with geometric $T$, $\kappa=L=1$ and the volumes of Section 3. [paper_figures/plots.py](paper_figures/plots.py) reads the saved CSVs; Figures 1 and 4 are schematics, with the example polynomials printed in Figure 4. The data schema, parameters and hashes are checked before drawing. These commands reproduce the published slices, not a general three-dimensional scan. See the [figure-to-data index](../data/paper_v3/README.md#asset-index).

Compare regenerated CSV values after decompression and SVG drawing content independently of timestamps and internal IDs. PNG pixel identity additionally depends on fonts and the rendering environment. The supplied manuscript images and data are retained as reference assets.

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

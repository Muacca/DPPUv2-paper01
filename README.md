# DPPUv2 Paper 01: Topology-Dependent Phase Classification of Effective Potentials in Einstein-Cartan Gravity

[日本語版 README](README_ja.md)

![Phase diagram example](data/image/ALL/phase_matrix_0010.png "Phase diagram example")

## Paper Content

- [Official PDF (English, built from LaTeX)](DPPUv2-paper01.pdf)
- [Draft Japanese version (Markdown)](20260101-01_paper_01_v2_2.md)

## Directory Structure

- `LaTeX/` - LaTeX manuscript and compiled PDF
  - `main.tex` - Main file
  - `sections/` - Section TeX files
  - `appendices/` - Appendix TeX files
  - `figures/` - Figures
- `data/` - Data for the paper
- `script/` - Data processing and visualization scripts
  - For script details, see [script/README](script/README.md)

## Building the LaTeX Document

### Basic Build Commands

Navigate to the LaTeX directory and run pdflatex.
To properly resolve cross-references (`\ref`, `\label`, etc.) in LaTeX, **you need to compile at least twice**:

- **1st pass**: Writes label information to `.aux` file
- **2nd pass**: Reads references from `.aux` file and resolves them

If references appear as `??`, compile once more.

```powershell
cd LaTeX
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

### Output Files

Upon successful build, the following files are generated:

- `main.pdf` - Final PDF file (approximately 3.5 MB, 79 pages)
- `main.aux` - Auxiliary file (cross-reference information)
- `main.log` - Compilation log
- `main.out` - Hyperref outline information

### Troubleshooting

#### If Errors Occur

Check errors:
```powershell
# Review log file
cat main.log | Select-String -Pattern "Error|Warning" | Select-Object -Last 20
```

#### Clean Build

Remove all auxiliary files before rebuilding:
```powershell
cd LaTeX
Remove-Item *.aux, *.log, *.out, *.synctex.gz, *.fdb_latexmk, *.fls -ErrorAction SilentlyContinue
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

## Author

Muacca

## License

See LICENSE file in the repository root.


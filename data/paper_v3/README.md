# Figure and table data

The paper uses geometric torsion $T$. These files use geometric $V,\eta$; the Python input dictionary is $V_\tau=-V,\eta_\tau=-\eta$. Scale $s$ means $r$ on the sphere and isotropic torus, or $R$ on Nil. The actual sphere radius is $a=r/2$. All numerical datasets use $\kappa=L=1$. Table 1 and Section 3 define the volumes and Nil lattice.

| Compressed CSV | Content |
|---|---|
| variant_curves.csv.gz | 1,201 scales per geometry, FULL/TT/positive REE and polynomial terms |
| representative_curves.csv.gz | 601 scales per geometry and each $\eta=-3,0,3$ |
| phase_grid.csv.gz | 115,515 geometry/variant/parameter rows; geometric axes and analytic labels |
| scaling.csv.gz | 121 values of $\lvert V\rvert$ for each valid-minimum example |
| representative_values.csv.gz | Table 4 extrema and heights |
| local_invariants.csv.gz | Table 5 local invariants |
| parameters.json | Parameters, grids, volumes, formula IDs, source-code and data hashes |

CSV files are UTF-8 with header rows, compressed with gzip. They can be read with Python's standard-library gzip and csv modules. Empty extrema or height fields mean that the strict minimum does not exist, not that its radius or height is zero. Phase labels I/II/III refer to Section 4. Type I height is a barrier; Type II height is the depth from the origin limit.

These figure datasets have their own descriptive columns; they are not the eight-column scanner output described in Appendix C. The phase slices are not a full default three-dimensional sweep. FULL agrees at $\theta=0,1,2$; TT and positive REE share the diagnostic dataset.

To reproduce a value with the attached scripts, pass the negated pair $(V_\tau,\eta_\tau)$ to the potential functions in script/DPPUv2_parameter_scan_v3.py. Tables 2 and 3 give polynomial coefficients and root conditions. For a curve, evaluate $P(As^3+Bs^2+Cs)$ at the saved scales. For a phase, classify positive stationary roots using $D=B^2-3AC$, treating $A=0$ separately. Table 5 gives local invariants before volume integration.

## Asset index

The same images appear in the Japanese Markdown and English PDF. PNG files use 360 dpi; SVG companions preserve vector curves. Phase colors are blue (I), orange (II), and gray (III). Analytic boundaries are overlaid on the 0.1-spaced phase cells.

| Figure | Asset | Dataset |
|---|---|---|
| 1 | [fig01_pipeline](../../LaTeX/figures/v3/fig01_pipeline.png) | Polynomial/classifier schematic |
| 2 | [fig02_variants](../../LaTeX/figures/v3/fig02_variants.png) | variant_curves.csv.gz |
| 3 | [fig03_decomposition](../../LaTeX/figures/v3/fig03_decomposition.png) | variant_curves.csv.gz |
| 4 | [fig04_extrema](../../LaTeX/figures/v3/fig04_extrema.png) | Polynomials printed in panels |
| 5 | [fig05_full_phases](../../LaTeX/figures/v3/fig05_full_phases.png) | phase_grid.csv.gz |
| 6 | [fig06_diagnostics_04_05](../../LaTeX/figures/v3/fig06_diagnostics_04_05.png) | phase_grid.csv.gz |
| 7 | [fig07_diagnostics_08_09](../../LaTeX/figures/v3/fig07_diagnostics_08_09.png) | phase_grid.csv.gz |
| 8 | [fig08_representatives](../../LaTeX/figures/v3/fig08_representatives.png) | representative_curves.csv.gz |
| 9 | [fig09_scaling](../../LaTeX/figures/v3/fig09_scaling.png) | scaling.csv.gz |

Tables 1–3 and 6–7 summarize symbolic definitions and criteria; Tables 4–5 use the datasets above. These are off-shell results along stated paths; action heights are not physical masses or decay rates.

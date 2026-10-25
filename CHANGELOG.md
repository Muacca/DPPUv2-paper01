# Revision history

## v3 — revision of v2 (2026-10-25)

This revision corrects the Nieh–Yan (NY) density and global volume normalization, and replaces the associated potential, phase-classification and vacuum interpretations. The date is the manuscript revision date.

### Action and torsion conventions

- Corrected the curvature-side NY density from $\frac14\varepsilon^{abcd}R_{abcd}$ to $\frac12\varepsilon^{abcd}R_{abcd}$. The torsion-side coefficient remains $1/4$. On the static homogeneous ansatz studied in the paper, the two contributions coincide, so the complete NY density vanishes and the FULL reduced potential is independent of $\theta$.
- Adopted geometric torsion $T$ in the main text. The appendices define the code-input array $\tau=-T$ and its parameter map $V_\tau=-V$, $\eta_\tau=-\eta$, including the reflected input domains. The code's contortion and connection formulas and input API are retained. This change of notation and convention does not itself change the potential.
- Distinguished the complete NY action from the TT and positive-REE diagnostic actions. The two diagnostic densities coincide on the stated ansatz; their effects must not be attributed to the complete action.

### Global geometry and normalization

The common $S^1$ coordinate has period 1, and $L$ denotes its physical circumference.

| Background | Adopted geometry | Four-volume |
|---|---|---|
| $S^3\times S^1$ | Round radius $a=r/2$ | $\pi^2 Lr^3/4$ |
| $T^3\times S^1$ | Spatial angles of period $2\pi$; circumferences $2\pi R_i$ | $(2\pi)^3 LR_1R_2R_3$ |
| $\mathrm{Nil}^3\times S^1$ | Heisenberg lattice $(2\pi,2\pi,4\pi^2)$ with twisted face identifications | $(2\pi)^4 LR^3$ |

At the same inputs, the volume correction multiplies the $S^3$ action, potential and derivatives by $1/8$, and the $T^3$ quantities by $1/(2\pi)$. The Nil volume is retained with an explicit compact quotient. These volume changes preserve radius extrema and reduced phase labels; the NY coefficient correction above changes the potentials and their classification.

### Conclusions and interpretation

- Withdrew the FULL action's NY-coupling-dependent phase boundaries and stabilization claims. With torsion coefficients and $L$ fixed, FULL has a strict $S^3$ radius minimum without a barrier for $V\ne0$ and $\lvert\eta\rvert<2$. It has no strict positive-radius minimum along isotropic or fixed-shape $T^3$ expansion, or uniform Nil expansion. Varying one torus radius while holding the other two fixed is treated as a distinct variation path.
- Retained barrier-separated minima only under the stated conditions for the diagnostic actions. Replaced the unsupported FULL-minimum and scaling examples accordingly.
- Reframed the classifications as off-shell extrema along specified variation paths. A radius minimum does not establish a full Einstein–Cartan vacuum, dynamical stability, or topology-only vacuum selection. Stated separately the obstruction to a smooth, nondegenerate full vacuum on closed real Euclidean $S^3\times S^1$ for source-free standard EC gravity without a cosmological constant and with constant NY coupling.
- Removed direct identifications of integrated-action heights with physical energy, mass or decay rates. Such an interpretation requires additional dynamical and normalization data.

### Manuscripts, computation and figures

- Revised the Japanese Markdown and English LaTeX/PDF consistently, including the title, action definitions, torsion dictionary, variation paths and interpretation limits.
- Renamed the Japanese Markdown chapters and appendices to stable filenames without date or version suffixes, and updated their links.
- Updated the calculation engine and its consumers for the corrected coefficients and volumes. Replaced finite-window searches for supported cubic potentials with analytic stationary-root classification; removed two-decimal input rounding and explicitly handled degenerate or numerically unresolved boundaries.
- Replaced the paper's illustrations and numerical summaries with nine figures and seven tables based on the revised equations. Supplied their input domains, parameters and numerical datasets in [data/paper_v3](data/paper_v3/README.md). Figure axes and datasets use geometric-torsion parameters. The phase datasets are specified slices, not a full default three-dimensional sweep.
- Updated the acknowledgments to disclose computational and plotting-code assistance used for this revision.
- Added an optional [independent calculation notebook](supplement/nieh_yan/README.md) deriving the local NY identity, its component coefficients and the T/τ dictionary directly from Cartan's equations, with English and Japanese explanations.

### Earlier results

Pre-revision images in [LaTeX/figures](LaTeX/figures/), historical datasets and saved notebook outputs are preserved but are not evidence for the v3 conclusions. Some generating inputs and versions remain unidentified; earlier scan statistics and version-comparison claims are not revalidated by this revision. Existing checkpoints and saved outputs are not automatically converted to the corrected coefficients or volume conventions.

The correction is tracked in [Issue #4](https://github.com/Muacca/DPPUv2-paper01/issues/4).

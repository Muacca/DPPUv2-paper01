# Nieh–Yan identity and the T/τ dictionary

[日本語](README_ja.md) · [Open the calculation notebook](nieh_yan_identity.ipynb)

This optional companion to Paper 01 derives the local Nieh–Yan expressions and the relation between geometric torsion $T$ and the code input $\tau$. The paper contains the definitions, claims and arguments in Section 2 and Appendices A–B; the notebook provides an executable way to follow that algebra.

## Calculation

The notebook contains English and Japanese explanations and uses exact symbolic arithmetic. It constructs the three coframes for $S^3$, $T^3$ and $\mathrm{Nil}^3$, solves Cartan's equation for a metric-compatible connection, and computes torsion and curvature as exterior forms. It then evaluates both sides of

$$
d(e^a\wedge T_a)=T^a\wedge T_a-e^a\wedge e^b\wedge R_{ab}.
$$

The component expansion gives the factors $1/4$ for TT and $1/2$ for REE. For the static homogeneous ansatz, all three geometries yield

$$
n_{\mathrm{TT}}=n_{\mathrm{REE}}=-\frac{4V\eta}{s},
\qquad n_{\mathrm{FULL}}=0.
$$

Here $s=r$ for the sphere, $`s=R_1`$ for the torus ansatz and $s=R$ for Nil; the sphere's physical radius is $r/2$. REE denotes the positive curvature term; FULL is TT minus REE. The notebook also compares the directly solved connection with both contortion formulas and establishes

$$
\tau=-T,\qquad V_\tau=-V,\qquad \eta_\tau=-\eta.
$$

Input intervals must be reflected along with the parameters. The notebook explains the endpoint ordering and gives an exact mixed-torsion example, as well as the axial-only and vector-only limits.

## Run

Use Python 3 with SymPy and Jupyter Notebook. From the repository root:

```powershell
python -m pip install sympy notebook
jupyter notebook supplement/nieh_yan/nieh_yan_identity.ipynb
```

Choose **Restart Kernel and Run All Cells**. The distributed notebook has no saved cell outputs. The density matrix has rows S3, T3 and Nil3 and columns TT, REE, FULL and $`d(e^a\wedge T_a)/\mathrm{vol}_4`$. Empty dictionaries in the residual displays denote zero exterior forms; the component-comparison matrix is zero.

The calculation imports only SymPy and Python's standard library. It reads no data files and does not call the paper's calculation engines. It can also be run after copying the notebook to a separate directory.

## Relation to the paper

- [English paper](../../DPPUv2-paper01.pdf): Section 2 and Appendices A–B give the action, connection conventions, component coefficients and torsion dictionary.
- [Japanese Section 2](../../section_02_ja.md), [Appendix A](../../appendix_A_ja.md) and [Appendix B](../../appendix_B_ja.md): the corresponding definitions and derivations.
- [Section 3](../../section_03_ja.md) and [conventions](../../script/docs/CONVENTIONS.md): coframe scales, integration domains and global volumes.

The calculation holds the parameters constant along the frame. It concerns local differential forms; the pointwise vanishing of FULL uses this particular ansatz, whereas the NY identity is general. Global normalization is specified separately in the paper. This notebook does not perform phase scans or solve the full field equations.

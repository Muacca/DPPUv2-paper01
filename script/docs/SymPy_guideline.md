# SymPy Implementation Guidelines — DPPUv2 Engine

- **Target:** AI Coding Assistants & Contributors
- **Context:** Einstein–Cartan action with a constant Nieh–Yan coupling on static homogeneous $S^3\times S^1$, $T^3\times S^1$ and compact $\mathrm{Nil}^3\times S^1$

This document gives the implementation rules for the symbolic engine `DPPUv2_engine_core_v3.py` and the three topology runners in this directory. Index order, frame, torsion and global-geometry conventions are fixed in [CONVENTIONS](CONVENTIONS.md); the paper states the same definitions in Section 2 and Appendix A.

⇒ [日本語版](SymPy_guideline_ja.md) | [Geometric Conventions](CONVENTIONS.md)

-----

## 1. Optimization Rules (Engineering)

Iron rules to prevent "computation time explosion" (taking over 2 hours) in SymPy integration calculations, ensuring completion within seconds.

The current runners use left-invariant frames with constant components and multiply the homogeneous density by the global volume (step `E4.12`). The rules below apply whenever explicit integration or large intermediate expressions are introduced.

### Rule 1.1: Strict use of `expand()` + `cancel()` Strategy

Do not use high-cost functions like `simplify()` on the integrand immediately before performing integration (`integrate`). Instead, perform "expansion and cancellation".

  * **Don't:**
    ```python
    density = simplify(density)  # NG: PROHIBITED. Factorizing huge expressions has extreme computational cost.
    result = integrate(density, x)
    ```
  * **Do:**
    ```python
    density = cancel(expand(density))  # OK: RECOMMENDED. Converting to a sum of polynomials induces term-by-term integration.
    result = integrate(density, x)
    ```

### Rule 1.2: Suppress Intermediate Simplification

When performing multiple integrations (e.g., $\phi$ integration $\to$ $\theta$ integration), do not apply excessive `simplify` to intermediate results. Limiting it to `cancel` and re-applying `expand` just before the final integration is faster.

-----

## 2. Implementation Rules for Tensor Components

### Rule 2.1: Use Frame Components in the Stored Index Order

The engine works in an orthonormal frame with the Euclidean metric $`g_{ab}=\delta_{ab}`$ (`metric_frame`, default `Matrix.eye(dim)`). Raising and lowering frame indices is then the identity, and components may be used directly. Array positions must follow [CONVENTIONS Section 2](CONVENTIONS.md#2-indices-and-array-index-order):

* `C[a,b,c]` $`=C^a{}_{bc}`$
* `Gamma[a,b,c]` $`=\Gamma^a{}_{bc}`$, where the last index is the direction of differentiation
* `Riemann[a,b,c,d]` $`=R^a{}_{bcd}`$

For a non-orthonormal basis or another metric, lower all indices with the metric, permute them, and raise them again. A bare array transposition is not an index operation.

### Rule 2.2: Distinguish the Input Array $\tau$ from Geometric Torsion $T$

`torsion_tensor` stores the input array $\tau$. The torsion obtained from the Cartan equation is the geometric torsion $T=-\tau$, stored as `geometric_torsion_tensor`. The contortion step `E4.5` applies the code formula to $\tau$:

```python
# Step E4.5: contortion of the input array tau
for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            val = (tau[a, b, c] + tau[b, c, a] - tau[c, a, b]) * Rational(1, 2)
```

In geometric torsion the same contortion is $`K_{abc}[T]=-\frac12(T_{abc}+T_{bca}-T_{cab})`$ (paper Eq. M2-03). Convert the array and both parameters together, $`\tau=-T`$, $`\eta_\tau=-\eta`$, $`V_\tau=-V`$. Do not reverse $K$ alone or only one torsion parameter.

### Rule 2.3: Torsion Consistency Check Including Structure Coefficients

In a noncoordinate frame, the torsion components contain the structure coefficients:

$$
T^a{}_{bc}=C^a{}_{bc}+\Gamma^a{}_{cb}-\Gamma^a{}_{bc}.
$$

After constructing the EC connection, confirm that this reproduces $-\tau$ with a mismatch count of **0**:

```python
dim = engine.data['dim']
C = engine.data['structure_constants']
Gamma_EC = engine.data['connection_EC']
tau = engine.data['torsion_tensor']
mismatch = sum(
    1
    for a in range(dim) for b in range(dim) for c in range(dim)
    if simplify(C[a, b, c] + Gamma_EC[a, c, b] - Gamma_EC[a, b, c] + tau[a, b, c]) != 0
)
assert mismatch == 0
```

The difference `Gamma_EC[a,b,c] - Gamma_EC[a,c,b]` alone omits $C$. It is not a valid check on $S^3$ or Nil, where $C\ne0$.

### Rule 2.4: Built-in Checks

* Step `E4.3` verifies metric compatibility $`\Gamma_{abc}+\Gamma_{bac}=0`$ for the Levi–Civita connection and raises an exception on violation.
* Step `E4.7` verifies $`R_{abcd}=-R_{bacd}`$ and $`R_{abcd}=-R_{abdc}`$ with the verdicts `PROVED_ZERO`, `WITNESS_NONZERO` and `UNPROVED`. An unproved residual is not treated as zero.
* Pair exchange $`R_{abcd}=R_{cdab}`$ is not required for a connection with torsion (paper Eq. B-02).

-----

## 3. Conventions Used by the Engine

### 3.1 Torsion and Curvature

Torsion and curvature are defined by two-forms, with $`\omega^a{}_b=\Gamma^a{}_{bc}e^c`$ ([CONVENTIONS Section 6](CONVENTIONS.md#6-torsion-and-curvature)):

$$
\begin{aligned}
T^a&=de^a+\omega^a{}_b\wedge e^b=\frac12T^a{}_{bc}\thinspace e^b\wedge e^c,\\
R^a{}_b&=d\omega^a{}_b+\omega^a{}_c\wedge\omega^c{}_b=\frac12R^a{}_{bcd}\thinspace e^c\wedge e^d.
\end{aligned}
$$

### 3.2 Connection

$$
\begin{aligned}
\Gamma^{\mathrm{EC}}_{abc}&=\Gamma^{\mathrm{LC}}_{abc}+K_{abc},\qquad
\Gamma^{\mathrm{LC}}_{abc}=\frac12\left(C_{abc}+C_{cba}-C_{bac}\right),\\
K_{abc}[T]&=-\frac12\left(T_{abc}+T_{bca}-T_{cab}\right)
=\frac12\left(\tau_{abc}+\tau_{bca}-\tau_{cab}\right).
\end{aligned}
$$

The Levi–Civita part is the general Koszul formula of [CONVENTIONS Section 5](CONVENTIONS.md#5-levi-civita-connection-general-koszul-implementation-in-dppuv2-engine). It does not assume a bi-invariant metric.

### 3.3 Curvature Components

With frame components constant in each runner, step `E4.7` uses

$$
R^a{}_{bcd}=\Gamma^a{}_{ec}\Gamma^e{}_{bd}-\Gamma^a{}_{ed}\Gamma^e{}_{bc}+\Gamma^a{}_{be}C^e{}_{cd}.
$$

The last term comes from $de^c$ in a noncoordinate frame and must not be dropped (paper Appendix G).

-----

## 4. Torsion Ansatz and Mode Decomposition Rules

### 4.1 Mode Definitions

| Mode | Input components |
|---|---|
| `Mode.AX` | Axial part only ( $`V_\tau=0`$ ) |
| `Mode.VT` | Vector-trace part only ( $`\eta_\tau=0`$ ) |
| `Mode.MX` | Both parts |

### 4.2 Components

Step `E4.4` builds the input array from the runner parameters `eta` $`=\eta_\tau`$, `V` $`=V_\tau`$ and the scale `r`. The scale alias is $r$ on $S^3$, $`R_1`$ on $T^3$, and $R$ on Nil; call it $s$.

$$
\begin{aligned}
\tau^{\mathrm{ax}}_{ijk}&=\frac{2\eta_\tau}{s}\varepsilon_{ijk}\quad(i,j,k\in\lbrace 0,1,2\rbrace ),\\
\tau^{\mathrm{vec}}_{abc}&=\frac13\left(\delta_{ac}(V_\tau)_b-\delta_{ab}(V_\tau)_c\right),\qquad
(V_\tau)_a=(0,0,0,V_\tau).
\end{aligned}
$$

The geometric torsion $T=-\tau$ has the same form with $`\eta=-\eta_\tau`$ and $`V=-V_\tau`$ (paper Eq. M2-05). Runners supply the metric, structure constants, volume and parameters in steps `E4.1`–`E4.2`; they do not construct the torsion array themselves.

-----

## 5. Nieh–Yan Densities and Action Variants

### 5.1 Densities

With $`\varepsilon_{0123}=+1`$ in the Euclidean frame, step `E4.10` computes

$$
N_{\mathrm{TT}}=\frac14\varepsilon^{abcd}T^e{}_{ab}T_{ecd},\qquad
N_{\mathrm{REE}}=\frac12\varepsilon^{abcd}R_{abcd},\qquad
N_{\mathrm{FULL}}=N_{\mathrm{TT}}-N_{\mathrm{REE}}.
$$

The TT factor $1/4$ comes from the two torsion two-forms and the REE factor $1/2$ from the single curvature two-form. $`N_{\mathrm{TT}}`$ is quadratic in torsion, so evaluating it with $\tau$ gives the same value. On the ansatz of the three runners,

$$
N_{\mathrm{TT}}=N_{\mathrm{REE}}=-\frac{4V\eta}{s}=-\frac{4V_\tau\eta_\tau}{s},\qquad N_{\mathrm{FULL}}=0.
$$

### 5.2 Variant Selection

| `NyVariant` | NY term in the Lagrangian | Role |
|---|---|---|
| `NyVariant.FULL` | $`\theta N_{\mathrm{FULL}}`$ | Complete NY action |
| `NyVariant.TT` | $`\theta N_{\mathrm{TT}}`$ | Diagnostic action |
| `NyVariant.REE` | $`+\theta N_{\mathrm{REE}}`$ | Diagnostic action with a positive sign; not the negative component of FULL |

On this ansatz TT and REE give the same potential. Their barriers and $\theta$ dependence are not results of the complete NY action.

### 5.3 Implementation

`step_E4_10_nieh_yan_frame` stores `ny_density_TT`, `ny_density_Ree` and `ny_density_full`, and selects the adopted density by the `ny_variant` argument of `BaseFrameEngine`.

-----

## 6. Action and Reduced Potential

### 6.1 Action

$$
L=\frac{\mathcal R_{\mathrm{EC}}}{2\kappa^2}+\theta N,\qquad
S=L\times\mathrm{Vol},\qquad
U=-S.
$$

Steps `E4.11`, `E4.12` and `E4.13` construct $L$, $S$ and $U$. The volume `total_volume` is defined by each runner ([CONVENTIONS Section 9](CONVENTIONS.md#9-global-geometry-and-integration-domains)). $U$ is stored as `potential`; the log labels it "Effective Potential V(r)".

### 6.2 Obtaining the Potential

```python
from DPPUv2_engine_core_v3 import Mode, NyVariant
from DPPUv2_runner_S3S1_v3 import S3S1Engine

engine = S3S1Engine(Mode.MX, NyVariant.FULL)
engine.run()
U = engine.get_effective_potential_function()
# Arguments: (scale, V_tau, eta_tau, theta_NY, L, kappa)
```

* The scale argument is `r` on $S^3$, `R1` on $T^3$ and `R` on Nil.
* The torsion arguments are code inputs. For the paper's geometric values, pass $`V_\tau=-V`$ and $`\eta_\tau=-\eta`$.
* On $T^3$ the stored potential keeps `R2` and `R3`. Substitute the variation path into `engine.data['potential']` (for example $`R_2=\alpha R_1`$, $`R_3=\beta R_1`$) before calling `get_effective_potential_function()`.
* `get_potential_decomposition()` returns the coefficients of powers of the scale.

-----

## 7. Classification of Radius Extrema

Step `E4.14` records derivatives only. Use the analytic classifier in `DPPUv2_parameter_scan_v3.py`:

* `get_cubic_coefficients(topology, ny_variant, V, eta, theta)` returns $(a,b,c)$ of $U=ar^3+br^2+cr$ for the registered isotropic $T^3$, sphere and uniform Nil paths, with code inputs.
* `analyze_cubic_stability(a, b, c)` classifies the positive-radius domain (paper Section 4, Table 3).

| Type | Definition on $r>0$ |
|---|---|
| Type-I | A positive-radius local maximum followed by a strict local minimum |
| Type-II | A strict local minimum without a preceding positive-radius maximum |
| Type-III | No strict local minimum: monotone, flat, runaway or stationary-inflection cases |

* The Type-I height is the action difference from the maximum to the minimum. The Type-II depth is measured from the origin limit $U(0)=0$. Do not plot them as the same kind of barrier.
* If $`b^2-3ac\le64\epsilon_{\mathrm{float}}\max(b^2,\lvert 3ac\rvert)`$, no strict minimum is resolved. This does not prove that the discriminant is zero; analytic inequalities take precedence.
* Callables without cubic coefficients fall back to a bounded numerical search, limited by its window and resolution.

These are off-shell radius extrema along fixed variation paths. They do not establish a full Einstein–Cartan vacuum, dynamical stability, nucleation, physical masses or decay rates (paper Sections 5 and 7).

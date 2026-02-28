# SymPy Implementation Guidelines — dppu Engine

- **Target:** AI Coding Assistants & Contributors
- **Context:** Einstein-Cartan Theory on Curved Spacetime (S³×S¹, T³×S¹, Nil³×S¹)

This document establishes the implementation guidelines for symbolic computation (SymPy) and theoretical physics conventions for the `dppu` package. The following rules must be strictly adhered to.

⇒ [日本語版](SymPy_guideline_ja.md) | [Geometric Conventions](CONVENTIONS.md)

-----

## 1. Optimization Rules (Engineering)

Iron rules to prevent "computation time explosion" (taking over 2 hours) in SymPy integration calculations, ensuring completion within seconds.

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

## 2. Theoretical Implementation Rules

Iron rules for maintaining consistency in tensor operations within Curved Spacetime.

### Rule 2.1: Prohibition of Direct Index Manipulation (Robust Method)

In environments with non-diagonal metrics or where $g_{\mu\nu} \neq 1$, swapping array indices directly (e.g., `T[mu, nu, lam]`) is not equivalent to physical tensor index manipulation (raising/lowering).

  * **Don't:**
    ```python
    # NG: PROHIBITED. Risk of referencing physically incorrect components.
    term = T_tensor[mu, nu, lam]
    ```

  * **Do:** Always manipulate indices via the metric tensor $g_{\mu\nu}$.
    1.  Lower all indices to create the fully covariant form $T_{\lambda\mu\nu}$.
    2.  Perform index permutation.
    3.  Raise indices using the metric as necessary.

#### Optimization for Orthonormal Frame Basis

In an **Orthonormal Frame Basis**, since the metric is the identity matrix ( $g_{ab} = \eta_{ab} = \text{diag}(1,1,1,1)$ or $\text{diag}(-1,1,1,1)$ ), the computational overhead of raising and lowering indices should be omitted. Direct component calculation should be used for speed.

However, the sign pattern $(+1, +1, -1)$ defined by the physical definition (Hehl 1976) must be strictly observed.

```python
# ============================================================
# Golden Logic for Contortion (Frame Basis / DPPUv2 Standard)
# ============================================================
# Assumption: Metric is diagonal/identity (Orthonormal Frame)
# Therefore T^a_bc and T_abc behave identically in code logic.

K_tensor = MutableDenseNDimArray.zeros(dim, dim, dim)

for a in range(dim):
    for b in range(dim):
        for c in range(dim):
            # Formula: K_abc = (1/2)(T_abc + T_bca - T_cab)
            # Note: Using T[a,b,c] directly as T_abc
            
            term = (T_tensor[a, b, c] + T_tensor[b, c, a] - T_tensor[c, a, b])
            
            val = term * Rational(1, 2)
            
            if val != 0:
                K_tensor[a, b, c] = cancel(expand(val))
```

### Rule 2.2: Consistency Check

After constructing the EC connection ($\Gamma_{\text{EC}}$), the following verification code must be executed to confirm that the mismatch count is **0**.

```python
# Torsion Consistency Check
T_verify = Gamma_EC[lam, mu, nu] - Gamma_EC[lam, nu, mu] # Hehl definition
mismatch = count(simplify(T_verify - T_original) != 0)
assert mismatch == 0
```

-----

## 3. Standard Conventions (Theoretical)

To prevent confusion during paper writing, adhere to the **Hehl (1976) Standard**.

### 3.1 Torsion Definition

Definition of the Torsion Tensor $T^\lambda_{\ \mu\nu}$:
$$T^\lambda_{\ \mu\nu} \equiv \Gamma^\lambda_{\ \mu\nu} - \Gamma^\lambda_{\ \nu\mu}$$
(Note: In the frame basis of this engine, torsion components are extracted from the torsion 2-form $T^a = de^a + \omega^{a}{}\_b\wedge e^b$ defined in Section 6 of CONVENTIONS, via the coefficient comparison $T^a = \frac{1}{2}T^{a}{}\_{bc}\,e^b\wedge e^c$.)

### 3.2 Contortion Formula

The Verified Formula for Contortion $K^\lambda_{\ \mu\nu}$ consistent with the above Torsion definition:

$$K_{\lambda\mu\nu} = \frac{1}{2} \left( T_{\lambda\mu\nu} + T_{\mu\nu\lambda} - T_{\nu\lambda\mu} \right)$$

  * Sign Pattern: **$(+1, +1, -1)$**
  * Note: Apply this formula to $T_{\lambda\mu\nu}$ (where all indices are lowered).

### 3.3 Einstein-Cartan Connection

$${\Gamma_{\text{EC}}}^\lambda_{\ \mu\nu} = {\Gamma_{\text{LC}}}^\lambda_{\ \mu\nu} + K^\lambda_{\ \mu\nu}$$

  * $\Gamma_{\text{LC}}$: Levi-Civita Connection (Christoffel symbols)
  * $K$: Contortion

-----

## 4. Torsion Ansatz and Mode Decomposition Rules

The torsion tensor on the $M^3 \times S^1$ minisuperspace ansatz is specified via three modes.

### 4.1 Mode Definitions

| Mode | Physical component | Parameters |
|---|---|---|
| `Mode.AX` | Axial component (T1) only | $\eta \neq 0$, $V = 0$ |
| `Mode.VT` | Vector-trace component (T2) only | $\eta = 0$, $V \neq 0$ |
| `Mode.MX` | Both T1 and T2 | $\eta \neq 0$, $V \neq 0$ |

### 4.2 Physical Correspondence

- **T1 (Axial component):** Dual to axial vector $S^\mu = (\eta/r)(0,0,0,1)$. For spatial indices $a,b,c \in \{0,1,2\}$: $T_{abc} = (2\eta/r)\,\varepsilon_{abc}$.
- **T2 (Vector-trace component):** Dual to vector $V_\mu = V\,\delta^3_\mu$ ($\tau$-component only). $T_{abc} = \frac{1}{3}(\delta_{ac}V_b - \delta_{ab}V_c)$.

### 4.3 Implementation Rule

Use `construct_torsion_tensor(mode, r, eta, V, metric, dim)` from `dppu/torsion/ansatz.py`.
Do **not** construct $T_{abc}$ by hand.

-----

## 5. Nieh-Yan Topological Term Variants

### 5.1 Nieh-Yan Decomposition

The full Nieh-Yan density:
$N = N_{\mathrm{TT}} - N_{\mathrm{Ree}},$
$N_{\mathrm{TT}} = \frac{1}{4}\varepsilon^{abcd}T^{e}{}\_{ab}T_{ecd},\qquad N_{\mathrm{Ree}} = \frac{1}{4}\varepsilon^{abcd}R_{abcd}.$

### 5.2 Variant Selection

| `NyVariant` | Density used |
|---|---|
| `NyVariant.TT` | $N_{\mathrm{TT}}$ only |
| `NyVariant.REE` | $N_{\mathrm{Ree}}$ only |
| `NyVariant.FULL` | $N_{\mathrm{TT}} - N_{\mathrm{Ree}}$ (canonical) |

### 5.3 Implementation

All three variants are computed in pipeline step `E4.10`. Select the variant via the `ny_variant` argument in the engine's `__init__`. See `dppu/torsion/nieh_yan.py`.

-----

## 6. Extended Lagrangian and Weyl Coupling Constant $\alpha$

### 6.1 Action Form

$$S = \int L \times \mathrm{Vol},\qquad
L = \frac{R}{2\kappa^2} + \theta_{\mathrm{NY}}\times N + \alpha\times C^2.$$

| Parameter | Meaning |
|---|---|
| $\kappa$ | Einstein-Cartan gravitational coupling |
| $\theta_{\mathrm{NY}}$ | Nieh-Yan coupling (topological) |
| $\alpha$ | Weyl coupling (conformal invariant term) |

For $\alpha \leq 0$, Theorem 1 guarantees protection of the stable vacuum. For $\alpha > 0$, Theorem 2 guarantees $V_{\rm eff} \to -\infty$.

### 6.2 Obtaining the Effective Potential

```python
engine.run()
V_func = engine.get_effective_potential_function()
# Call signature:
# V_func(r, V_param, eta, theta_NY, L, kappa, epsilon, alpha) -> float
```

$V_{\rm eff} = -S$. Extracted in pipeline step `E4.13`.

### 6.3 Implementation

See `compute_lagrangian()` in `dppu/action/lagrangian.py` and `dppu/action/potential.py`.

-----

## 7. Numerical Optimisation Strategy (Phase Atlas Search)

### 7.1 Two-Stage Strategy

**Stage 1: Brute-force grid search**

`scipy.optimize.brute` with `Ns` points per axis over the $(r, \varepsilon)$ 2D grid, to locate the basin of the global minimum.

**Stage 2: Multi-start L-BFGS-B refinement**

Top- $N$ candidates from Stage 1 are used as starting points for `scipy.optimize.minimize` (L-BFGS-B, `ftol=1e-8`) to obtain a high-precision minimum.

### 7.2 Stability Classification

After locating the minimum at $(r^\*, \varepsilon^\*)$ , classify as follows:

| Class | Condition | Physical meaning |
|---|---|---|
| Type-I | $V(r)$ rises as $r \to 0$ (barrier present) | Stable vacuum (nucleation barrier exists) |
| Type-II | $V(r)$ decreases monotonically as $r \to 0$ (no barrier) | Spontaneous nucleation possible |
| Type-III | No local minimum in the physical region | Unstable configuration |

Use `analyze_stability()` from `dppu/action/stability.py`.

### 7.3 Important Notes

- For $\alpha > 0$, the optimiser attaches to the search boundary $(r \to 0^+,\,\varepsilon \to -1^+)$ (`converged = False`). This is not an optimisation failure; it correctly reports that the potential is unbounded below within the search range.
- For $\mathrm{Nil}^3$ flat-limit confirmation, the search range is extended to $\varepsilon_{\rm max} = 5.0$, wider than the default (see `scripts/paper02/`).

-----

## 8. Self-Duality (SD) Diagnostic Rules

### 8.1 Hodge Dual of Curvature

$$(*R)^{ab}{}_{cd} = \frac{1}{2}\varepsilon_{cdef}\,R^{ab,ef},$$

where $\varepsilon_{cdef}$ is the Levi-Civita symbol in the frame basis (use `levi_civita_4d()` from `dppu/utils/levi_civita.py`).

### 8.2 Pontryagin Inner Product and SD Residuals

$$E_{RR} = \langle R, R\rangle = R_{abcd}R^{abcd},\qquad
P = \langle R, *R\rangle = R_{abcd}(*R)^{abcd}.$$

| Condition | Physical state |
|---|---|
| SD residual $< \varepsilon_{\rm SD}$ and $\|R\| > \varepsilon_R$ | Self-dual |
| ASD residual $< \varepsilon_{\rm SD}$ and $\|R\| > \varepsilon_R$ | Anti-self-dual |
| $P = 0$ (Proposition 1) | Chiral equilibrium (algebraic identity on $M^3 \times S^1$) |

### 8.3 Usage

```python
from dppu.curvature.self_duality import SDExtensionMixin
SDExtensionMixin.attach_to(engine)          # Dynamically attach methods
R = engine.get_R_ab_cd_numerical(params_dict)
diag = engine.evaluate_sd_status(params_dict)
# Check diag['P_RstarR'] (expected value: 0)
```

See `dppu/curvature/self_duality.py` and `dppu/curvature/pontryagin.py`.



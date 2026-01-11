# DPPUv2 Implementation Guidelines (Engineering & Theoretical Guidelines)

- **Version:** 1.0
- **Target:** AI Coding Assistants
- **Context:** Einstein-Cartan Theory on Curved Spacetime (S3xS1, etc.)

This document establishes the implementation guidelines for symbolic computation (SymPy) and theoretical physics conventions for the DPPUv2 project. The following rules must be strictly adhered to.

⇒ [日本語版](DPPUv2_SymPy_guideline_v1_ja.md)

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

In an **Orthonormal Frame Basis**, since the metric is the identity matrix ($g_{ab} = \eta_{ab} = \text{diag}(1,1,1,1)$ or $\text{diag}(-1,1,1,1)$), the computational overhead of raising and lowering indices should be omitted. Direct component calculation should be used for speed.

However, the sign pattern $(+1, +1, -1)$ defined by the physical definition (Hehl 1976) must be strictly observed.

```python
# ============================================================
# Golden Logic for Contortion (Frame Basis / v3 Standard)
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
                K_tensor[a, b, c] = simplify(val)
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
(Note: When calculating from a tetrad, $T^a_{\mu\nu} = \partial_\nu e^a_\mu - \partial_\mu e^a_\nu$ corresponds to this).

### 3.2 Contortion Formula

The Verified Formula for Contortion $K^\lambda_{\ \mu\nu}$ consistent with the above Torsion definition:

$$K_{\lambda\mu\nu} = \frac{1}{2} \left( T_{\lambda\mu\nu} + T_{\mu\nu\lambda} - T_{\nu\lambda\mu} \right)$$

  * Sign Pattern: **$(+1, +1, -1)$**
  * Note: Apply this formula to $T_{\lambda\mu\nu}$ (where all indices are lowered).

### 3.3 Einstein-Cartan Connection

$${\Gamma_{\text{EC}}}^\lambda_{\ \mu\nu} = {\Gamma_{\text{LC}}}^\lambda_{\ \mu\nu} + K^\lambda_{\ \mu\nu}$$

  * $\Gamma_{\text{LC}}$: Levi-Civita Connection (Christoffel symbols)
  * $K$: Contortion



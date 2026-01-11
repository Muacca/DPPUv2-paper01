# CONVENTIONS (DPPUv2 v3)

This document fixes the **geometric, index, and sign conventions** shared between `BaseFrameEngine (v3)` and each topology runner.
**All runners must define `metric_frame` and `structure_constants` according to these conventions.**

⇒ [日本語版](DPPUv2_engine_core_v3_CONVENTIONS_ja.md)

## 1. Scope and Assumptions

* All quantities treated here are expressed as components on a **Frame (Orthonormal Basis)**.
* `metric_frame` is the frame metric $(g_{ab})$, which defaults to $(g_{ab}=\delta_{ab})$ (`Matrix.eye(dim)`).
* The curvature component calculation in the v3 engine currently assumes a situation where **frame directional derivatives are unnecessary**.
  Specifically, the runner must adopt a setting where structure constants $(C^a{}_{bc})$ and connection coefficients $(\Gamma^a{}_{bc})$ are treated as "constants with respect to the frame" (e.g., left-invariant frames).
  * **Note:** This design is **rational and efficient** for handling homogeneous spaces like $S^3 \times S^1$ (Lie groups) or Nil manifolds. However, note that this assumption may become a bottleneck if handling general curved spacetimes with low symmetry in the future.

## 2. Indices and Array Index Order

Array storage is fixed as follows:

* Structure constants: `C[a,b,c] = C^a_{bc}`
* Connection coefficients: `Gamma[a,b,c] = Γ^a_{bc}`
* Riemann curvature: `Riemann[a,b,c,d] = R^a_{bcd}`

Meaning of indices:

* $(a)$: Output (superscript) component
* $(b,c,d)$: Input (subscript) components
  In particular, $(\Gamma^a{}_{bc})$ corresponds to $(\nabla_{E_c} E_b = \Gamma^a{}_{bc} E_a)$ (the last index $(c)$ is the "direction of differentiation").

## 3. Definitions of Frame, Coframe, and Structure Constants (Most Important)

Let the dual of the frame be $(\{E_a\})$ and the coframe (1-forms) be $(\{e^a\})$.

### 3.1 Structure Equations of Coframe (Fixed)

$$
de^a = \frac12 C^a{}_{bc} e^b\wedge e^c,
\qquad C^a{}_{bc} = - C^a{}_{cb}.
$$

### 3.2 Commutation Relations of Dual Frame (Equivalent)

The above definition is equivalent to:
$$
[E_b, E_c] = - C^a{}_{bc} E_a.
$$

> Note: Many textbooks adopt $([E_b,E_c]=+f^a{}_{bc}E_a)$.
> This project adopts the convention **$(C^a{}_{bc}=-f^a{}_{bc})$ relative to that $(f^a{}_{bc})$**.

### 3.3 Runner Implementation Rules (Recommended)

* **Do not manually input C**. If possible, specify `de^a` on the runner side, extract `C^a_{bc}` by coefficient comparison, and put it into `self.data['structure_constants']=C`.
* At a minimum, automatically check that $C^a_{bc}$ is **antisymmetric in b,c**.

## 4. Connection (Spin Connection) and Metric Compatibility

Define the connection 1-form as:
$$
\omega^a{}_b = \Gamma^a{}_{bc} e^c
$$

Metric compatibility (Lorentz connection / Orthogonal connection) is fixed as a specification:
$$
\omega_{ab} = -\omega_{ba}
\quad(\Leftrightarrow\quad
\Gamma_{abc} = -\Gamma_{bac})
$$
where $(\Gamma_{abc} = g_{ad}\Gamma^d{}_{bc})$.

## 5. Levi-Civita Connection (General Koszul Implementation in v3 Engine)

When the frame is orthonormal and the above structure constant conventions are adopted, the Levi-Civita connection is calculated in the v3 engine using the following **General Koszul Formula**:

$$
\Gamma^a{}_{bc}
= \frac12\Big(
C^a{}_{bc} + C^c{}_{ba} - C^b{}_{ac}
\Big).
$$

(This is the form when the sign of the commutation relation in Section 3.2 is adopted.)

**Important Notes:**

1. This formula **does not assume a bi-invariant metric**.
   It functions correctly as the Levi-Civita connection on left-invariant frames, even for non-bi-invariant cases like Nil³.

2. In special cases where structure constants are totally antisymmetric in lower indices $C_{abc} = -C_{bac} = -C_{acb}$ (like SU(2)),
   this formula reduces to $\Gamma^a_{bc} = \frac{1}{2} C^a_{bc}$.

3. The engine automatically verifies **metric compatibility** $\Gamma_{abc} + \Gamma_{bac} = 0$ after calculation.
   If this is violated, it immediately throws an exception as an implementation error.

## 6. Torsion and Curvature

Torsion 2-form:
$$
T^a = de^a + \omega^a{}_b \wedge e^b,
\qquad
T^a = \frac12 T^a{}_{bc} e^b\wedge e^c.
$$

Curvature 2-form:
$$
R^a{}_b = d\omega^a{}_b + \omega^a{}_c\wedge \omega^c{}_b,
\qquad
R^a{}_b = \frac12 R^a{}_{bcd} e^c\wedge e^d.
$$

## 7. Curvature Component Formula (Form used by Engine)

The v3 engine currently uses the following form for $R^a_{bcd}$:

$$
R^a{}_{bcd}
=
\Gamma^a{}_{ec}\Gamma^e{}_{bd}
-\Gamma^a{}_{ed}\Gamma^e{}_{bc}
-\Gamma^a{}_{be} C^e{}_{cd}.
$$

> Important: Generally, a frame directional derivative term
> $(E_c(\Gamma^a{}_{bd}) - E_d(\Gamma^a{}_{bc}))$
> appears here, but the v3 engine does not explicitly handle it.
> Therefore, the runner must adopt a setting where **$(\Gamma)$ is treated as constant in the frame direction**, such as with left-invariant frames.

## 8. Mandatory Self-Checks (Consistency the Runner Must Satisfy)

The runner must satisfy the following (failure implies definition inconsistency with the engine):

1. Antisymmetry of Structure Constants:
   $$
   C^a{}_{bc} + C^a{}_{cb} = 0.
   $$

2. Metric Compatibility (Orthogonal Connection):
   $$
   \omega_{ab} + \omega_{ba} = 0.
   $$

3. Antisymmetry of Riemann (Target of engine's strict check):
   $$
   R_{ab cd} = -R_{ba cd},\qquad
   R_{ab cd} = -R_{ab dc}.
   $$



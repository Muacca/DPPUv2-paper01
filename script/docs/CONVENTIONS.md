# CONVENTIONS — dppu Engine Geometric Conventions

This document fixes the **geometric, index, and sign conventions** shared between `BaseFrameEngine` and each topology runner in the `dppu` package.
**All runners must define `metric_frame` and `structure_constants` according to these conventions.**

⇒ [日本語版](CONVENTIONS_ja.md) | [SymPy ガイドライン](SymPy_guideline.md)

## 1. Scope and Assumptions

* All quantities treated here are expressed as components on a **Frame (Orthonormal Basis)**.
* `metric_frame` is the frame metric $`(g_{ab})`$, which defaults to $`(g_{ab}=\delta_{ab})`$ (`Matrix.eye(dim)`).
* The curvature component calculation in the DPPUv2 engine currently assumes a situation where **frame directional derivatives are unnecessary**.
  Specifically, the runner must adopt a setting where structure constants $`(C^{a}{}_{bc})`$ and connection coefficients $`(\Gamma^{a}{}_{bc})`$ are treated as "constants with respect to the frame" (e.g., left-invariant frames).
  * **Note:** This design is **rational and efficient** for handling homogeneous spaces like $S^3 \times S^1$ (Lie groups) or Nil manifolds. However, note that this assumption may become a bottleneck if handling general curved spacetimes with low symmetry in the future.

## 2. Indices and Array Index Order

Array storage is fixed as follows:

* Structure constants: `C[a,b,c] = C^a_{bc}`
* Connection coefficients: `Gamma[a,b,c] = Γ^a_{bc}`
* Riemann curvature: `Riemann[a,b,c,d] = R^a_{bcd}`

Meaning of indices:

* $(a)$: Output (superscript) component
* $(b,c,d)$: Input (subscript) components
  In particular, $`(\Gamma^{a}{}_{bc})`$ corresponds to $`(\nabla_{E_c} E_b = \Gamma^{a}{}_{bc} E_a)`$ (the last index $(c)$ is the "direction of differentiation").

## 3. Definitions of Frame, Coframe, and Structure Constants (Most Important)

Let the dual of the frame be $`(\{E_a\})`$ and the coframe (1-forms) be $`(\{e^a\})`$.

### 3.1 Structure Equations of Coframe (Fixed)

$$
de^a = \frac12 C^a{}_{bc} e^b\wedge e^c,
\qquad C^a{}_{bc} = - C^a{}_{cb}.
$$

### 3.2 Commutation Relations of Dual Frame (Equivalent)

The above definition is equivalent to:

$$
[E_b, E_c] = - C^{a}{}_{bc} E_a.
$$

> Note: Many textbooks adopt $`([E_b,E_c]=+f^{a}{}_{bc}E_a)`$.
> This project adopts the convention **$`(C^{a}{}_{bc}=-f^{a}{}_{bc})`$ relative to that $`(f^{a}{}_{bc})`$**.

### 3.3 Runner Implementation Rules (Recommended)

* **Do not manually input C**. If possible, specify `de^a` on the runner side, extract `C^a_{bc}` by coefficient comparison, and put it into `self.data['structure_constants']=C`.
* At a minimum, automatically check that $`C^a_{bc}`$ is **antisymmetric in b,c**.

## 4. Connection (Spin Connection) and Metric Compatibility

Define the connection 1-form as:

$$
\omega^{a}{}_b = \Gamma^{a}{}_{bc} e^c
$$

Metric compatibility (Lorentz connection / Orthogonal connection) is fixed as a specification:

$$
\omega_{ab} = -\omega_{ba}
\quad(\Leftrightarrow\quad
\Gamma_{abc} = -\Gamma_{bac})
$$

where $`(\Gamma_{abc} = g_{ad}\Gamma^d{}_{bc})`$.

## 5. Levi-Civita Connection (General Koszul Implementation in DPPUv2 Engine)

When the frame is orthonormal and the above structure constant conventions are adopted, the Levi-Civita connection is calculated in the DPPUv2 engine using the following **General Koszul Formula**:

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

2. In special cases where structure constants are totally antisymmetric in lower indices $`C_{abc} = -C_{bac} = -C_{acb}`$ (like SU(2)),
   this formula reduces to $`\Gamma^a_{bc} = \frac{1}{2} C^a_{bc}`$.

3. The engine automatically verifies **metric compatibility** $`\Gamma_{abc} + \Gamma_{bac} = 0`$ after calculation.
   If this is violated, it immediately throws an exception as an implementation error.

## 6. Torsion and Curvature

Torsion 2-form:

$$
T^a = de^a + \omega^{a}{}_b \wedge e^b,
\qquad
T^a = \frac12 T^{a}{}_{bc} e^b\wedge e^c.
$$

Curvature 2-form:

$$
R^{a}{}_b = d\omega^a{}_b + \omega^a{}_c\wedge \omega^c{}_b,
\qquad
R^{a}{}_b = \frac12 R^{a}{}_{bcd} e^c\wedge e^d.
$$

### 6.1 Input array and geometric torsion

`torsion_tensor` stores the input parameter array tau. The geometric components are
`T_geom[a,b,c] = C[a,b,c] + Gamma[a,c,b] - Gamma[a,b,c] = -tau[a,b,c]`.
`geometric_torsion_tensor` exposes these components. The contortion formula is
`K[a,b,c]=(tau[a,b,c]+tau[b,c,a]-tau[c,a,b])/2`.

### 6.2 Nieh–Yan density

For Euclidean frame metric and epsilon(0,1,2,3)=+1, the component contractions are
`N_TT=epsilon*T*T/4`, `N_REE=epsilon*R/2`, and `N_FULL=N_TT-N_REE`.
REE denotes the diagnostic action term `+theta*N_REE`.
For the static homogeneous ansatz of the three runners,
`TT=REE=-4*V*eta/r` (using each runner's scale alias), hence `FULL=0`.

## 7. Curvature Component Formula (Form used by Engine)

The DPPUv2 engine currently uses the following form for $`R^a_{bcd}`$:

$$
R^{a}{}_{bcd} = \Gamma^{a}{}_{ec}\Gamma^{e}{}_{bd}-\Gamma^{a}{}_{ed}\Gamma^{e}{}_{bc}+\Gamma^{a}{}_{be} C^{e}{}_{cd}.
$$

> Important: Generally, a frame directional derivative term
> $`(E_c(\Gamma^{a}{}_{bd}) - E_d(\Gamma^{a}{}_{bc}))`$
> appears here, but the DPPUv2 engine does not explicitly handle it.
> Therefore, the runner must adopt a setting where **$(\Gamma)$ is treated as constant in the frame direction**, such as with left-invariant frames.

## 8. Mandatory Self-Checks (Consistency the Runner Must Satisfy)

The runner must satisfy the following (failure implies definition inconsistency with the engine):

1. Antisymmetry of Structure Constants:

$$
C^{a}{}_{bc} + C^{a}{}_{cb} = 0.
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

---

## 9. Weyl Tensor and Conformal Scalar

### 9.1 Weyl Tensor Definition (4D)

$$
C_{abcd} = R_{abcd} - \frac{1}{2}(g_{ac}R_{bd} - g_{ad}R_{bc} - g_{bc}R_{ad} + g_{bd}R_{ac}) + \frac{R}{6}(g_{ac}g_{bd} - g_{ad}g_{bc}).
$$

Key properties:
- **Traceless**: $`C^a{}_{bad} = 0`$ (holds for all index pairs)
- **Conformally invariant**: invariant under $`g_{ab} \to \Omega^2 g_{ab}`$
- **Conformal flatness criterion**: $`C_{abcd} = 0 \Leftrightarrow`$ conformally flat

Note on frame basis (orthonormal): since $`g_{ab} = \delta_{ab}`$, raising and lowering indices is the identity operation; components $`C_{abcd}`$ and $C^{abcd}$ can be treated identically.

### 9.2 Weyl Scalar

$$
C^2 = C_{abcd}\thinspace C^{abcd} = \sum_{a,b,c,d} C_{abcd}^2.
$$

In the frame basis, the sum is computed directly without index raising.
$C^2 = 0$ holds for isotropic $S^3 \times S^1$ at $\varepsilon = 0$ and for $T^3 \times S^1$ identically.

### 9.3 Engine Implementation

- Module: `dppu/curvature/weyl.py`
  - `compute_weyl_tensor(R_abcd, Ricci, R_scalar, metric, dim)` → $`C_{abcd}`$
  - `compute_weyl_scalar(C_abcd, metric_inv, dim)` → $C^2$
- Pipeline step: `E4.3b` (immediately after the Levi-Civita curvature computation)

---

## 10. Squashed Homogeneous Spaces and the $\varepsilon$-Parameter

### 10.1 Definition of Squashing

A volume-preserving anisotropy deformation parameter $\varepsilon$ is introduced.
$\varepsilon = 0$ is the isotropic reference point. Physical range: $\varepsilon \in (-1, +\infty)$ ($\varepsilon = -1$ is a singularity where structure constants diverge).

### 10.2 Topology-Dependent Structure Constant Scaling

The base structure constants $`C^a{}_{bc}(\varepsilon=0)`$ on the left-invariant frame are scaled as follows.

**$S^3 \times S^1$ (SU(2)):**

| Frame index $a$ | Scale factor $`\lambda_a(\varepsilon)`$ |
|---|---|
| $`a \in \{0, 1\}`$ | $(1+\varepsilon)^{2/3}$ |
| $a = 2$ | $(1+\varepsilon)^{-4/3}$ |

Volume preservation: $`\lambda_0\lambda_1\lambda_2 = (1+\varepsilon)^{2/3+2/3-4/3} = 1`$.

**$\mathrm{Nil}^3 \times S^1$ (Heisenberg group):**

| Non-trivial index $a$ | Scale factor |
|---|---|
| $a = 2$ (non-commutative component) | $(1+\varepsilon)^{-4/3}$ |

As $\varepsilon \to +\infty$, the structure constants vanish and the geometry asymptotes to flat $T^3$-like behaviour.

**$T^3 \times S^1$ (Abelian group):**

Structure constants are identically zero ($`C^a{}_{bc} = 0`$). The $\varepsilon$-deformation is not defined. $C^2 = 0$ holds identically, so the Weyl term vanishes (null test).

### 10.3 Physical Limits

| Limit | Physical meaning |
|---|---|
| $\varepsilon = 0$ | Isotropic $S^3$: $C^2 = 0$, stable vacuum of Paper I |
| $\varepsilon \to +\infty$ | $\mathrm{Nil}^3$ flat limit ($C^2 \to 0$ asymptotically) |
| $\varepsilon \to -1^+$ | Singularity of $S^3$ structure (physically inaccessible) |
| $\varepsilon = -2$ | Mathematical root of $C^2 = 0$, excluded since $\varepsilon < -1$ |

## 11. Global Geometry and Integration Domains

This section specifies the global backgrounds used by the three topology runners. The metric is the sum of the squares of the orthonormal coframe. All scales are positive. On the common circle, $u\in[0,1)$ and $`e^3=L\,du`$, so **L is the circumference**. The engine integrates a homogeneous scalar density by multiplying it by the Riemannian volume of the specified domain.

### 11.1 Round S³

The scale parameter r in $`C^i{}_{jk}=4\varepsilon_{ijk}/r`$ is twice the round-sphere radius $a=r/2$. An explicit SU(2) coframe is $e^i=-r\lambda^i/4$, where

$$
\begin{aligned}
\lambda^0&=\cos\psi\thinspace d\vartheta+\sin\psi\sin\vartheta\thinspace d\varphi,\\
\lambda^1&=-\sin\psi\thinspace d\vartheta+\cos\psi\sin\vartheta\thinspace d\varphi,\\
\lambda^2&=d\psi+\cos\vartheta\thinspace d\varphi.
\end{aligned}
$$

Use $0\le\vartheta\le\pi$, $0\le\varphi<2\pi$, $0\le\psi<4\pi$, with the usual polar degeneracies and coupled boundary identification $(\varphi+2\pi,\psi)\sim(\varphi,\psi+2\pi)$. These are not three independently periodic torus coordinates. The positive metric density is $r^3\sin\vartheta/64$; the displayed coframe has orientation opposite to $d\vartheta\wedge d\varphi\wedge d\psi$.

$$
\mathrm{Vol}(S^3)=\frac{\pi^2r^3}{4}=2\pi^2a^3,\qquad
\mathrm{Vol}(S^3\times S^1)=\frac{\pi^2Lr^3}{4},\qquad
R_{\mathrm{LC}}=\frac{24}{r^2}=\frac{6}{a^2}.
$$

### 11.2 T³

For $i=1,2,3$, take $x^i\in[0,2\pi)$ and $`e^{i-1}=R_i\,dx^i`$. The $`R_i`$ are coframe scale coefficients; the spatial circumferences are $`\ell_i=2\pi R_i`$.

$$
\mathrm{Vol}(T^3\times S^1)=(2\pi)^3LR_1R_2R_3
=L\ell_1\ell_2\ell_3.
$$

The scan uses $`R_1=R_2=R_3=r`$. A fixed-shape path $`R_2=\alpha r`$, $`R_3=\beta r`$ holds positive $\alpha,\beta$ fixed. Holding R2 and R3 themselves fixed defines a different variation path.

### 11.3 Compact Nil³

Use the Heisenberg group law and left-invariant coframe

$$
\begin{aligned}
(x,y,z)(x',y',z')&=(x+x',y+y',z+z'+xy'),\\
e^0&=R\thinspace dx,\qquad e^1=R\thinspace dy,\qquad e^2=R(dz-x\thinspace dy).
\end{aligned}
$$

The compact space is the left quotient $\Gamma\backslash H$ with

$$
\Gamma=\lbrace (2\pi j,2\pi k,4\pi^2l):j,k,l\in\mathbb Z\rbrace .
$$

Its generators identify coordinates by

$$
\begin{aligned}
(x,y,z)&\sim(x+2\pi,y,z+2\pi y),\\
(x,y,z)&\sim(x,y+2\pi,z),\\
(x,y,z)&\sim(x,y,z+4\pi^2).
\end{aligned}
$$

Each identification preserves the coframe. For periods $(a,b,c)=(2\pi,2\pi,4\pi^2)$, the horizontal-generator commutator is the central translation by $ab=c$, so $ab/c=1$. The coordinate box $[0,2\pi)\times[0,2\pi)\times[0,4\pi^2)$ is a fundamental domain with the twisted x-face gluing above. Integration uses the positive density $`R^3\,dx\,dy\,dz`$.

$$
de^2=-\frac1R e^0\wedge e^1,\qquad C^2{}_{01}=-\frac1R,\qquad
\mathrm{Vol}(\Gamma\backslash H\times S^1)=(2\pi)^4LR^3.
$$

The base-torus circumferences are $2\pi R$ and the central-fiber circumference is $4\pi^2R$. The lattice and its shape are held fixed as R varies. This explicit quotient realizes the local structure constants and the stated global volume together; its central coordinate does not have period $2\pi$.

### 11.4 Input Parameters and Physical Lengths

The runner and consumer APIs use the local scales r/R1/R2/R3/R and the input-array parameters $`\eta_\tau,V_\tau`$. The geometric torsion dictionary remains $`\eta_g=-\eta_\tau`$, $`V_g=-V_\tau`$.

If the same S³ axial torsion is written as $2\chi/a$ instead of $`2\eta_g/r`$, use $a=r/2$ and $`\chi=\eta_g/2`$. For T³ circumference variables, $`2\zeta/\ell_1=2\eta_g/R_1`$ requires $`\zeta=2\pi\eta_g`$. Parameter intervals, sampling grids and extrema transform with these dictionaries. V remains a signed coefficient and acquires no additional factor.

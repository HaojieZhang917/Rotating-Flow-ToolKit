# Compressible BEK linear-stability refactor and operator audit

Date: 2026-09-11

## Scope

This audit covers the compressible rotating-disk spatial stability model used
by the first three cells of
`/home/zhj/Rotating-Flow-ToolKit/compress/BEK/roughness/CRD_STA_CACU.ipynb`
and compares it with Turkyilmazoglu & Uygun (2006), *Direct Spatial Resonance
in the Compressible Boundary Layer on a Rotating Disk*.

The verified reference is deliberately narrow: the local, parallel,
compressible von Karman boundary layer at `Ro=-1`. It is a spatial linear
stability calculation, not a nonlinear stability calculation and not a
global finite-radius calculation.

## Clean calculation graph

The refactored program separates six operations that were interleaved in the
notebook:

1. physical parameters: `Tw, Mr, Ro, R, beta, omega, gamma, Pr`;
2. dependent parameters: `Ma=Mr/R`, `Co=2-Ro-Ro^2`;
3. velocity and thermal base-flow BVPs plus endpoint checks;
4. sampling of every base variable on the stability coordinate `y`;
5. assembly of `(A0 + alpha*A1 + alpha^2*A2)q=0`;
6. balancing, targeted solution, backward-error verification and output.

The shared interface is
`work/src/CompressibleBEKLinearStability.jl`; the executable interface is
`work/scripts/run_compressible_bek_lsa_20260911.jl`.

## Paper-to-code definitions

The paper uses perturbations proportional to

\[
  \widehat\phi(y)\exp[i(\alpha r+\beta\theta-\omega t)]
\]

and obtains a quadratic spatial eigenproblem in complex \(\alpha\). The
paper's local rotational Mach number is \(M_r\), while the inherited primitive
operator expects

\[
  Ma=\frac{M_r}{R}.
\]

The five retained unknowns are ordered as

\[
 (\widehat\rho,\widehat u,\widehat v,\widehat w,\widehat T).
\]

Pressure and Chapman-law property perturbations are eliminated by

\[
 \widehat p=\frac{T\widehat\rho+\rho\widehat T}{\gamma Ma^2},\qquad
 \widehat\mu=\widehat T,\qquad
 \widehat\lambda=-\frac23\widehat T,\qquad
 \widehat k=\frac{\widehat T}{Pr}.
\]

This is algebraically valid at finite `Ma`, but the factors
`1/(gamma*Ma^2)` make the raw matrices badly scaled as `Ma -> 0`. The new
solver balances the polynomial matrices and judges modes using a polynomial
backward error.

## Findings and corrections

### 1. Stability-coordinate error in the notebook

The notebook used `interp(...,"phy")`. Equations (10)--(15) of the paper are
written in the Dorodnitsyn--Howarth similarity coordinate \(y\), so the
operator and differentiation matrices require `interp(...,"sim")`.

The physical-coordinate reconstruction is post-processing. Inserting it
between the base solution and the published stability operator would require
all derivative coefficients to be transformed as well.

### 2. Thermal BVP endpoint conditions were not imposed

The inherited `f_q` routine passed a multipoint boundary callback but used
`u[begin]` and `u[end]` rather than evaluating the solution at the interval
endpoints. A direct diagnostic returned

\[
  f(40)=0.021989\ldots,\qquad q(40)=0.085944\ldots
\]

instead of zero. The new implementation solves two genuine Julia
`TwoPointBVProblem`s and rejects the base state unless

\[
 \max(|f(0)|,|f(40)|,|q(0)-1|,|q(40)|)<10^{-8}.
\]

This defect is weak in the present `Mr=0.1, Tw=1` benchmark because the
viscous-heating correction is proportional to `Mr^2`, but it is unacceptable
for heated or higher-Mach calculations.

### 3. Spatial operator rows

For the strict von Karman slice, the compact tensor assembly contains the
five expected equations: continuity; radial, azimuthal and wall-normal
momentum; and energy. It retains pressure closure, Coriolis/curvature
couplings, normal viscous terms, heat conduction, base-property variation and
the Chapman-law disturbance terms.

At each validation grid, the compact tensor assembly was compared with the
independently expanded `Spatial_mode_BEK1` implementation after applying the
same boundary elimination. Relative discrepancies are at roundoff
(approximately `1e-16` for `A0`, `1e-24` for `A1`, and zero for `A2`). This
shows that the two code paths implement the same algebra; it is not by itself
an independent derivation of the general-Ro model.

### 4. Boundary conditions and order counting

The density-like variable is not fixed at the wall. The eliminated degrees
of freedom impose no-slip velocity and isothermal-temperature perturbations
at the wall and decay at the far boundary, plus the far density condition.
Nine eliminated degrees of freedom give the correct order count for the
five-field isothermal problem. The current boundary eliminator does not
implement an insulated thermal condition.

### 5. Artificial regularisation

The old compact constructor defaulted to an `O(1e-8)` density term in the
quadratic continuity block. The new interface uses exactly zero by default.
Nonzero regularisation is exposed only as a sensitivity diagnostic and is
recorded in the output metadata.

### 6. Status of the Ro extension

The 2006 paper directly verifies only the von Karman equations. The inherited
BEK extension uses `Co=2-Ro-Ro^2`, sign-reverses some base profiles depending
on `Ro`, and maps the special input `Ro=-1` to `+1` inside the stability
tensors. This convention reproduces the paper slice but is not a transparent
derivation of arbitrary-Ro compressible perturbation equations.

The new API therefore rejects `Ro != -1` unless the caller explicitly sets
`allow_unverified_ro=true`; such output is labelled
`unverified_bek_extension`. General-Ro output is not suitable for a manuscript
claim until the rotating-frame primitive equations are nondimensionalised and
linearised row by row and an independent benchmark is supplied.

## Verification case

The prescribed case is

\[
 M_r=0.1,\quad T_w=1,\quad Ro=-1,\quad R=285.36,\quad
 \beta=0.07759,\quad\omega=0,
\]

so `Ma=3.50434538828e-4`. The published `alpha_r=0.38482` is used as the
eigensolver target, not imposed as an answer. The final table is written to
`work/results/compressible_bek_lsa_validation_20260911/grid_convergence.csv`.

The corrected calculation gives:

| domain | N | alpha | backward error |
|---|---:|---:|---:|
| rational | 59 | `0.3851895715 + 0.0003923048im` | `1.53e-14` |
| rational | 69 | `0.3851895712 + 0.0003923072im` | `1.51e-14` |
| rational | 79 | `0.3851895692 + 0.0003923045im` | `7.54e-15` |
| finite, ymax=40 | 79 | `0.3851895722 + 0.0003922959im` | `3.79e-14` |
| finite, ymax=40 | 99 | `0.3851895702 + 0.0003923059im` | `3.16e-14` |

Across the three rational grids, the full spread is approximately `2.3e-9`
in the real part and `2.7e-9` in the imaginary part. The finite and rational
grids agree to about `1e-8` in the imaginary part. The thermal-BVP endpoint
error is `4.19e-17`; compact/expanded operator errors remain at roundoff.

The converged real part differs from the rounded literature target by about
`9.60e-4` relative, while the imaginary part is small (`3.92e-4`) but not
zero. Thus this fixed `(R,beta)` case is a well-resolved near-neutral
regression. It is not a recomputation of the exact neutral point: that requires
continuing one parameter and enforcing `imag(alpha)=0` while tracking the same
mode by eigenvector overlap.

Acceptance criteria are:

- thermal-BVP endpoint residual below `1e-8`;
- compact/expanded operator discrepancies near roundoff;
- QEP polynomial backward error below `1e-10`;
- eigenvalue stability under `N` and domain changes;
- an explicit distinction between a rounded near-neutral literature point
  and a newly recomputed exact neutral point.

The completed grid study gives

| grid | \(N\) | \(\alpha\) | polynomial backward error |
|---|---:|---:|---:|
| rational | 59 | \(0.3851895715+0.0003923048i\) | \(1.53\times10^{-14}\) |
| rational | 69 | \(0.3851895712+0.0003923072i\) | \(1.51\times10^{-14}\) |
| rational | 79 | \(0.3851895692+0.0003923045i\) | \(7.54\times10^{-15}\) |
| finite, \(y_{max}=40\) | 79 | \(0.3851895722+0.0003922959i\) | \(3.79\times10^{-14}\) |
| finite, \(y_{max}=40\) | 99 | \(0.3851895702+0.0003923059i\) | \(3.16\times10^{-14}\) |

Across the rational grids the spread is about \(2.3\times10^{-9}\) in the
real part and \(2.7\times10^{-9}\) in the imaginary part. The finite-domain
and rational results agree to approximately \(10^{-8}\). The recomputed mode
is therefore numerically converged for this formulation.

Its real part differs from the paper's rounded `0.38482` by approximately
`9.6e-4` relative, and the imaginary part is small but not zero. Thus this
parameter tuple is a well-reproduced *near-neutral* point. Determining the
discretely exact neutral point requires releasing one control parameter and
solving `imag(alpha)=0`; that continuation was intentionally outside this
validation-only calculation.

## Remaining work before arbitrary-Ro production runs

1. derive the general rotating-frame compressible BEK perturbation equations,
   including the definitions and signs of `Ro` and `Co`;
2. compare every continuous coefficient with the assembled matrix blocks;
3. construct a low-Mach continuation check that separates primitive-variable
   conditioning from physics;
4. validate at least one additional published compressible neutral point;
5. only then enable unrestricted `Ro` neutral-curve continuation.

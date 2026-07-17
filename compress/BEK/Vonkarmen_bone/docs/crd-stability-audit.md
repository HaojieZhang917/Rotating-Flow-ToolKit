# CRD_STA.jl matrix audit against Turkyilmazoglu and Uygun (2006)

## 1. Scope and conventions

This audit compares `CRD_STA.jl` with equations (10)-(15) of
Turkyilmazoglu and Uygun (2006), *Direct spatial resonance in the
compressible boundary layer on a rotating-disk*.

The strict paper-compatible parameter choices are

```text
Ro = -1
Co = 2
Ma = Mr / R
mu_B = T_B
lambda_B = -2 T_B / 3
k_B = T_B / sigma
rho_B = 1 / T_B
```

The paper uses the normal mode

```text
exp(i(alpha r + beta theta - omega t))
```

and the Dorodnitsyn-Howarth coordinate `y`. Therefore the stability
operator must use the `sim` base flow. The `phy` profile is only suitable
for physical-coordinate post-processing.

Define

```text
F = ubar_B, G = vbar_B, H = wbar_B,
U = alpha F + beta G - omega, D = d/dy.
```

The paper retains six primitive perturbations `(u,v,w,p,rho,T)`. The code
uses the state equation and Chapman law to eliminate

```text
p_hat = (T rho_hat + rho T_hat) / (gamma Ma^2),
mu_hat = T_hat,
lambda_hat = -2 T_hat / 3,
k_hat = T_hat / sigma.
```

The spatial code then uses the column ordering

```text
q_hat = (rho_hat, u_hat, v_hat, w_hat, T_hat)^T.
```

This reduction is algebraically equivalent to the paper for nonzero Mach
number, but it is not the same numerical formulation. In particular, it
introduces entries of size `R/(gamma*Ma^2)`, which are about `6.8e8` at
`R=440.88, Mr=0.3`.

## 2. Spatial matrix: equation-by-equation audit

### Continuity row: paper equation (10)

The code contains all five blocks:

```text
rho_hat: i R U + 2F + rho H' + rho H D
u_hat:   (i alpha R + 1) rho
v_hat:   i beta R rho
w_hat:   R rho (rho' + rho D)
T_hat:   0
```

Verdict: **consistent**. See `Spatial_mode_BEK` coefficient blocks
`Ta_14`, `A_14`, `B_14`, `C_14`, `C_13`, and `D_11`-`D_14`.

### Radial momentum row: paper equation (11)

After pressure and viscosity elimination, the required blocks are

```text
u_hat:
  i R rho U + rho F + (lambda+2T) alpha^2 + T beta^2
  + rho^2 H D - rho D^2

v_hat:
  -2 rho (G+1) + alpha beta (lambda+T)

w_hat:
  R rho^2 F' - i alpha [rho T' + (1+lambda rho)D]

rho_hat:
  rho F'' + i alpha R T/(gamma Ma^2)

T_hat:
  -rho D(rho F') - rho^2 F' D
  + i alpha R rho/(gamma Ma^2)
```

With `Ro=-1` converted internally to `Ro=1` while `Co=2` is retained,
`-rho*(2Ro*G+Co)` becomes `-2rho(G+1)`, exactly as in the paper.

Verdict: **consistent**.

### Azimuthal momentum row: paper equation (12)

The required blocks are

```text
u_hat:
  2 rho (G+1) + alpha beta (lambda+T)

v_hat:
  i R rho U + rho F + (lambda+2T) beta^2 + T alpha^2
  + rho^2 H D - rho D^2

w_hat:
  R rho^2 G' - i beta [rho T' + (1+lambda rho)D]

rho_hat:
  2F(G+1) + rho H G' + i beta R T/(gamma Ma^2)

T_hat:
  -rho D(rho G') - rho^2 G' D
  + i beta R rho/(gamma Ma^2)
```

Verdict: **consistent**.

### Axial momentum row: paper equation (13)

The required blocks are

```text
u_hat:
  -i alpha [rho lambda' + (1+lambda rho)D]

v_hat:
  -i beta [rho lambda' + (1+lambda rho)D]

w_hat:
  i R rho U + rho^2(H' + H D)
  - rho(2+lambda rho)D^2 + T(alpha^2+beta^2)

rho_hat:
  R rho [T' + T D]/(gamma Ma^2)

T_hat:
  -i rho(alpha F' + beta G')
  + R rho [rho' + rho D]/(gamma Ma^2)
```

Verdict for `Spatial_mode_BEK` and `Spatial_mode_BEK1`: **consistent**.

The temporal implementation has a separate error in this row; see section 4.

### Energy row: paper equation (14)

The code correctly includes:

```text
- velocity-gradient work in the u, v, and w blocks;
- base-temperature convection R rho^2 T' w_hat;
- pressure-work terms after state-equation elimination;
- heat conduction with k_B=T/sigma and k_hat=T_hat/sigma;
- viscous dissipation proportional to
  -(gamma-1) Ma^2 rho^2 (F'^2+G'^2) T_hat.
```

After combining the direct temperature convection and pressure work, the
coefficient of `i R rho U T_hat` becomes `1/gamma`, which is exactly what
the coefficient representation stores in `Ta_55`, `A_55`, and `B_55`.

Verdict: **consistent**.

## 3. Independent spatial implementation comparison

`Spatial_mode_BEK1` directly expands the polynomial

```text
L(alpha) = L0 + alpha L1 + alpha^2 L2.
```

`Spatial_mode_BEK` stores the coefficient tensors and `assemble_mat`
forms the same polynomial. An independent smooth-profile comparison gave

```text
relative difference L0 = 1.02e-12
relative difference L1 = 4.43e-24
relative difference L2 = 4.07e-9
```

The `L0` absolute differences occur in pressure-scaled blocks of order
`1e8`-`1e9` and are floating-point association differences. `L1` agrees
to roundoff. The only structural `L2` difference is the intentional
`1e-8` continuity-density regularization in `Vxx_14`.

Verdict: the two spatial implementations are **algebraically equivalent
apart from the documented regularization**.

## 4. `Timemode` correction status

The original temporal expansion contained three transcription errors: the
axial-momentum `beta*lambda'` multiplication operator and two energy-row
`omega` coefficients. Its boundary deletion also constrained wall density.

The active `Timemode` no longer maintains a second hand-written operator.
It is generated directly from the verified spatial coefficient tensors at
fixed `alpha`:

```text
B0 = L0(omega=0) + alpha L1 + alpha^2 L2,
B1 = i Ta,
B0 q = omega B1 q.
```

This also unifies the temporal and spatial variable ordering as
`(rho,u,v,w,T)`. The old formula is retained only as
`Timemode_legacy` for comparison and must not be used for new results.

The optimized regression test gives relative temporal matrix differences
below `1e-14` by construction.

## 5. Boundary conditions

### Spatial problem

The spatial problem has differential order nine after pressure
elimination. The code removes nine degrees of freedom:

```text
u=v=w=T=0 at y=0 and y=infinity,
rho=0 at y=infinity,
rho at the wall remains free.
```

This is consistent for an isothermal wall and does not prescribe wall
pressure. It also matches the required order count.

The implementation does not support an insulated wall. That case requires
`D T_hat=0` at the wall rather than deleting the wall temperature degree
of freedom.

### Temporal problem

The active temporal operator now calls the same nine-condition boundary
routine as the spatial operator. Wall density remains free and far-field
density is zero.

## 6. Difference from the paper's numerical method

The paper retains pressure and the state equation, and uses staggered
Chebyshev points: continuity on Gauss-Lobatto points and the remaining
equations on Gauss points. `CRD_STA.jl` eliminates pressure and uses one
collocated grid.

The formulations are algebraically equivalent for finite Mach number, but
they are not numerically identical. Pressure elimination produces very
large coefficients at low local Mach number and makes scaling important.

## 7. Type-II convergence diagnosis

The following test used

```text
R=440.88, beta=0.04672, omega=0,
Tw=1, Mr=0.3, Ma=Mr/R, sim coordinate.
```

Balanced generalized-Schur results with singular `A2` were

| N | alpha |
|---:|:---|
| 29 | 0.1271950322 - 0.0049904048i |
| 39 | 0.1326673575 + 0.0012648908i |
| 49 | 0.1326359919 + 0.0012450488i |
| 59 | 0.1326359364 + 0.0012451165i |
| 69 | 0.1326359369 + 0.0012451169i |

The physical target mode is therefore converged by about `N=49`-`69`.

When IAR requested four eigenpairs, it reported no convergence at
`N=39,49,59` because not all four Ritz pairs met `tol=1e-11`. The target
Type-II pair itself had a polynomial residual near `1e-13`. Requesting
several modes and treating the whole call as failed incorrectly labels the
target mode as non-converged.

The rational Chebyshev mapping becomes severely ill-conditioned:

| N | cond(interior D2) |
|---:|---:|
| 49 | 8.9e10 |
| 69 | 1.4e12 |
| 99 | 2.5e13 |
| 149 | 6.5e14 |
| 199 | 6.6e15 |

At `N=199`, double precision is near its useful limit. More points no
longer imply a more accurate eigenvalue.

The leading polynomial matrix has condition number about `1.39e8` when
the `1e-8` regularization is used. IAR also uses a random starting vector
unless `v` is supplied. At high `N`, changing the regularization or random
start produced eigenvalue drift of order `1e-6` even with small backward
residuals. This is an eigensolver-conditioning effect, not disappearance
of the Type-II mode.

## 8. Far-field truncation test

A finite-domain Chebyshev grid greatly improves differentiation-matrix
conditioning, but the Type-II far field is longer than the main velocity
profiles suggest.

At `N=99`:

| ymax | alpha |
|---:|:---|
| 20 | 0.1341809156 + 0.0020249290i |
| 30 | 0.1327347998 + 0.0012880151i |
| 40 | 0.1326431871 + 0.0012462581i |

Thus `ymax=20` is not sufficient for this Type-II mode. `ymax=40` is close
to the rational-domain result, while retaining much better conditioning.

## 9. Recommended numerical configuration

1. Use only the `sim` coordinate for the stability operator.
2. For the current rational mapping, use approximately `N=59`-`79`.
3. Alternatively use a finite domain near `ymax=40` and verify domain
   convergence before grid convergence.
4. Request one target eigenpair, or at most the two physical candidates.
5. Return and reuse the previous eigenvector as IAR's `v`; identify the
   mode by eigenvector overlap as well as eigenvalue distance.
6. Apply consistent left/right polynomial matrix equilibration.
7. IAR and generalized Schur can operate with singular `A2`; do not assume
   the `1e-8` term is mathematically required. If retained, expose it as a
   parameter and include a sensitivity test.
8. Use the active tensor-generated `Timemode`; do not use
   `Timemode_legacy` for results.
9. `Wcc_fun` now assigns both endpoint weights. Adjoint boundary conditions
   still require an independent physical validation before sensitivity work.

## 10. Optimized solver entry point

`solve_spatial_mode` now provides the recommended path. It defaults to a
singular leading matrix, deterministic initialization, one target mode,
consistent left/right polynomial equilibration, and returns the original
unscaled eigenvectors plus backward residuals.

For continuation, pass the previous returned eigenvector through
`initial_vector`. The notebook wrapper exposes this as `v0`.

The base-flow helper caches the `N`-independent velocity and thermal
profiles in each Julia session. Changing collocation resolution or domain
therefore no longer repeats the SciPy BVP and thermal auxiliary solves.

Notebook usage:

```julia
eigval,eigvec = eigsol(
    F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,omega,be,
    N_cheb,Ro,Co,D,D2,0.13,1,
)

# Continue the same mode at the next parameter point.
eigval_next,eigvec_next = eigsol(
    F,G,H,rho,lam,kappa,T,sigma,gamma,R_next,Ma_next,omega,be_next,
    N_cheb,Ro,Co,D,D2,eigval[1],1; v0=eigvec[:,1],
)

# Request residual and scaling diagnostics.
info = eigsol(
    F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,omega,be,
    N_cheb,Ro,Co,D,D2,0.13,1; return_info=true,
)
```

## 11. Material-property diagnostic switches

The active spatial and temporal operators expose two independent switches:

- `property_perturbations=false` removes only the Chapman-law terms
  `mu_hat=T_hat` and `kappa_hat=T_hat/sigma`. Thermodynamic temperature,
  density, pressure, and the energy equation remain active.
- `base_property_variation=false` freezes `mu_bar`, `lambda_bar`, and
  `kappa_bar` at their far-field values inside the perturbation operator.
  It does not recompute or replace the supplied base flow.

Recommended three-case comparison:

```julia
full = solve_spatial_mode(args...; target=alpha0)

no_property_perturbations = solve_spatial_mode(
    args...; target=alpha0,
    property_perturbations=false,
)

frozen_base_properties = solve_spatial_mode(
    args...; target=alpha0,
    property_perturbations=false,
    base_property_variation=false,
)
```

The first difference isolates material-property disturbances. The second
isolates the variable base transport coefficients while retaining the same
compressible base profiles. A model-consistent constant-property study would
additionally require recomputing the base flow with the same property law.

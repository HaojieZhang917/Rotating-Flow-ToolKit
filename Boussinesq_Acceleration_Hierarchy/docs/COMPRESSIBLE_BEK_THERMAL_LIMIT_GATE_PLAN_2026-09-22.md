# General-Ro compressible BEK thermal-limit Gate: model and study registration

Date: 2026-09-22. This records the model and validation sweep **before**
running the new calculation. The evidence is intended for solver validation,
not a manuscript-ready neutral curve.

## Parent scaling and model

Use signed `DeltaOmega=Omega_F-Omega_D`, `Ro=DeltaOmega/Omega`,
`Co=2-Ro-Ro^2`, `Mr=r*abs(DeltaOmega)/a_inf`, `Ma=Mr/R`, Chapman `mu=T`,
`k=T/Pr`, and Lingwood direct `(U,V,W)`. Normal distance is scaled with
`sqrt(nu/Omega)`. Radial and normal heat advection are of order
`DeltaOmega`, while normal heat conduction is of order `Omega`; their ratio
is `Ro`. Viscous heating normalized against conduction is proportional to
`Mr^2`, without an additional explicit `Ro`. With the temperature
decomposition `T=1-(gamma-1)Mr^2*f/2+(Tw-1)q`, the source equations are

```text
f'' - Pr*Ro*W*f' - 2Pr*Ro*U*f = 2Pr*(U'^2+V'^2)
q'' - Pr*Ro*W*q' = 0.
```

For exact `Ro=0` and nonzero thermal forcing, the equations become
`f''=2Pr*shear` and `q''=0`. The prescribed wall and ambient Dirichlet data
cannot both be attained on the half-line when shear is nonzero. The
far-field homogeneous exponent is `Pr*Ro*W_inf`; if it is nonnegative, the
current similarity thermal model has no decaying mode. The code shall reject
these cases by default and allow them only as `:finite_diagnostic` with an
explicit `evidence_scope`.

At fixed system-rotation Mach `Msystem=r*Omega/a_inf`, the matched path is
`Mr=abs(Ro)*Msystem`. Its exact `Ro=0` endpoint has `Mr=0`; the current
compressible pressure scaling `1/Ma^2` is then singular and must not be used
as an incompressible eigensolver. This is an explicit limit prescription,
not an exact compressible Ekman solution at finite `Mr`.

## Registered tests

- Test the default guard at `(Ro,Mr,Tw)=(0,0.1,1)`, `(0,0,1.1)`, and
  `(1,0.1,1)`, while retaining `:finite_diagnostic` as an explicit path.
- Test `Ro=-1` with `Mr=0.1,Tw=1` as an admissible half-line candidate;
  evaluate structural operator gates, directional QEP derivative, and a
  targeted mode's backward error.
- Compare the fixed-system-Mach path `Ro=-0.2,-0.1,-0.05` with
  `Msystem=0.3,Tw=1,R=116.3,beta=0.137`; record `Mr`, temperature amplitude,
  thermal far-field exponent, structural gates and mode roots where feasible.
- Because the thermal decay length grows as `1/abs(Pr*Ro*W_inf)`, run a
  follow-up source-profile domain check at `Ro=-0.2` with `L=40,80` and
  `Ro=-0.1,-0.05` with `L=40,80,160`. Compare `f` at fixed near-wall points,
  the physical temperature perturbation and `f(0.75L)/max|f|` before treating
  any near-zero-Ro point as numerically half-line converged.
- Following the first domain sweep, require at least ten far-field thermal
  e-folds by default and repeat the matched-system-Mach stability path with
  `L=80,160,320` for `Ro=-0.2,-0.1,-0.05`, respectively. The initial
  `L=40` modes remain diagnostics, not accepted half-line results.
- Re-evaluate an explicit `Ro=0,Mr=0.1,Tw=1` finite-domain diagnostic to
  ensure the previous numerical root remains reproducible but is not
  misclassified as a physical half-line result.
- At `Ro=1,Mr=0.1,Tw=1` run the operator structure and QEP algebra gates
  only in explicit finite-domain-diagnostic scope; its thermal outflow flag
  prevents the root from being promoted to a physical neutral result.
- Repeat the `Ro=-0.05,L=320` matched-Mach case at `N=59` to check a
  representative eigenvalue's grid sensitivity after thermal-domain repair.
- Make the public `run_compressible_bek_lsa_20260911.jl` CLI expose the new
  thermal model and decay-length Gate, and replace its call to the retired
  compact/expanded audit with the active structural and QEP-derivative Gates.
- On the accepted `Ro=-1,L=40` and `Ro=-0.1,L=160` thermal cases, substitute
  the sampled `(U,V,W,T,rho,f,q)` into the separately coded primitive
  continuity, two momentum, state and energy residuals over `y<=12`, using
  finite Chebyshev grids `N=120` and `N=240`, respectively. Report residuals
  as diagnostics even if they miss a convergence threshold.

Use `gamma=1.4`, `Pr=0.72`, rational Chebyshev `N=49` for the path and
`N=69` for the matched diagnostic, `ymax=40`; an optional `N=59` numerical
check may be added after the first pass. Output directory:
`work/results/compressible_bek_thermal_limit_gate_20260922/`. The preserved
legacy backend and all `baselines/` are read-only dependencies.

# Compressible BEK thermal BVP correction and test

Date: 2026-09-15

## Change

In `work/src/CompressibleBEKLinearStability.jl`, the thermal source solver was
changed from

```text
q'' + Pr*Ro*W*q' = 0
f'' + Pr*Ro*W*f' + 2Pr*Ro*U*f = 2Pr*(U'^2+V'^2)
```

to the equations obtained from the unified primitive energy equation,

```text
q'' - Pr*Ro*W*q' = 0
f'' - Pr*Ro*W*f' - 2Pr*Ro*U*f = 2Pr*(U'^2+V'^2).
```

The source solver's tridiagonal convention is `y'' - transport*y' -
reaction*y = rhs`, so this corresponds to `transport=Pr*Ro*W` and
`reaction=2Pr*Ro*U`. The preserved `RotatingDiskFlow/src/CRD_STA.jl` backend
and `baselines/` were not changed.

## Primitive residual test

The test uses `Ro=(-1,0,1)`, `Mr=0.3`, `Tw=1`, `Pr=0.72`, `gamma=1.4`,
`N=(80,120,160,200)`, and evaluates the pre-DH primitive equations on both
`y<=12` and `y<=39`. At `N=120`, the literature-compatible energy residuals
are:

| flow | core absolute residual | term-balance ratio | extended absolute residual |
|---|---:|---:|---:|
| von Karman, `Ro=-1` | `8.09e-8` | `1.78e-6` | `8.09e-8` |
| Ekman, `Ro=0` | `3.68e-7` | `2.59e-6` | `3.68e-7` |
| Bodewadt, `Ro=1` | `3.15e-7` | `2.98e-6` | `3.15e-7` |

The old-sign equations now show residuals of `1.87e-3` (von Karman) and
`6.84e-2` (Bodewadt), confirming that the reduction is testing an independent
equation rather than reproducing the source solver's own equation.

## Remaining far-field limitation

The local differential residuals pass, but a finite-domain solution still has
to be checked against the half-line boundary condition. For wall heating
`(Mr,Tw)=(0,0.8)`, corrected `q` gives:

- von Karman: `q(1)=0.678`, `q(10)=2.82e-3`, `q(20)=4.84e-6`, so the decay is
  physically located near the wall;
- Bodewadt: `q(1)=1`, `q(30)=0.99994`, `q(39)=0.622`, so the transition is
  still pinned to the artificial endpoint `y=40`;
- Ekman: `q=1-y/40`, reflecting the exact `Ro=0` degeneration `q''=0` and
  the absence of a decaying half-line solution with both Dirichlet values.

Thus the corrected thermal differential equation passes at all three
endpoints, while the Bodewadt wall-temperature problem and the exact Ekman
thermal limit require a separate half-line or nonsimilar treatment.

## (M_r=0.3) stability check

The corrected source was used in continuous fixed-beta spatial-mode searches.
The physically anchored Ekman branch at `beta=0.132` gives the `N=49`
candidate `R=119.6237`, `alpha=0.49804+0i`; re-solving at `N=69` at that same
`R` gives `alpha=0.49833+2.27e-4 i`. The branch is close to neutral and close
to the previous Ekman candidate, but still needs an independent `N=69`
neutral solve.

The Bodewadt branch is strongly grid-sensitive: the `N=49` search gives
`R=56.5816`, `alpha=0.43724+0i`, while evaluating that point at `N=59` and
`N=69` gives respectively `alpha=0.46241-0.09045i` and
`0.52317-0.01498i`. It therefore does not provide a converged compressible
neutral parameter. At the supplied incompressible reference tuple
`(R,beta)=(27.4,0.115)` the `Mr=0.01` mode is
`alpha=0.48170-0.00180i`, close to the incompressible reference
`alpha_r=0.487`; at `Mr=0.3` after continuation to `beta=0.127` it is
`alpha=0.54080-0.20592i`, clearly non-neutral.

## Decision

The heat-equation sign error is fixed and the independent primitive energy
residual Gate passes for the differential equations. The complete BEK
compressible base is not yet fully validated because Bodewadt wall-heating
and the exact `Ro=0` thermal limit remain far-field singular/domain-dependent,
and the Bodewadt `Mr=0.3` neutral branch is not grid converged. The current
neutral values remain diagnostic candidates until those two issues are
resolved.

Outputs:

- `work/results/compressible_bek_baseflow_primitive_residual_gate_v3_20260915/`
- `work/results/compressible_bek_lingwood_anchored_neutral_v3_20260915/`
- `work/results/compressible_bek_corrected_fixed_beta_neutral_v3_20260915/`

## Domain-length and (Ro=0) follow-up

`source_profiles` now accepts `ymax` and `prepare_case` passes the numerical
domain length through to the thermal source grid. A domain sweep at `N=80`
was run for `L_y=40,60,80`.

For Bodewadt wall heating, `q(L_y-1)=0.621518` at all three lengths, while
`q(10)` remains essentially 1. The thermal transition therefore translates
with the artificial endpoint and is not a converged half-line solution. For
the viscous-heating stability case, the fixed reference-point mode changes
from `0.50028-0.20810i` (`L_y=40`) to `0.49993-0.20949i` (`L_y=60`) and
`0.50528-0.20851i` (`L_y=80`).

Re-solving the Bodewadt neutral condition independently at each length gives
`R=56.58165` (`L_y=40`), `56.61959` (`L_y=60`), and `65.70576`
(`L_y=80`). The last value is on a different grid-sensitive branch: at its
fixed `R`, `N=59` and `N=69` give imaginary parts `-0.05744` and `-0.10136`.
The (M_r=0.3) Bodewadt neutral parameter is therefore unresolved.

At exact `Ro=0`, `prepare_case` now marks any nontrivial thermal input
(`Mr>0` or `Tw!=1`) with `evidence_scope=:ro0_thermal_degenerate`. The finite
interval solution remains available as an explicit diagnostic, but is not
silently labeled a regular derived BEK base state. The domain sweep confirms
the reason: `q(1)` changes from `0.975` to `0.9833` to `0.9875` as
`L_y=40,60,80`, while the exact equation is `q''=0`.

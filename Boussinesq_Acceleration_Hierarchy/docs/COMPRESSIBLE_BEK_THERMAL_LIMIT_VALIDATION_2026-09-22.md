# Compressible BEK thermal-limit implementation and validation

Date: 2026-09-22. Model and sweep were registered beforehand in
`docs/COMPRESSIBLE_BEK_THERMAL_LIMIT_GATE_PLAN_2026-09-22.md`.

## Decision from the parent scaling

The adopted BEK convention is `Ro=DeltaOmega/Omega`,
`Mr=r*abs(DeltaOmega)/a_inf`, `Co=2-Ro-Ro^2`, with normal length
`sqrt(nu/Omega)`. Mean-flow thermal advection is therefore multiplied by
`Ro`; dissipation is proportional to `Mr^2`, not an additional `Ro^2`.
The current `source_profiles` already implements

```text
f'' - Pr*Ro*W*f' - 2Pr*Ro*U*f = 2Pr*(U'^2+V'^2),
q'' - Pr*Ro*W*q' = 0.
```

No coefficient was changed simply to force a neutral root. The implemented
change is to enforce the domain of validity of this thermal similarity
system. At exact `Ro=0`, `q''=0` cannot satisfy `q(0)=1,q(infinity)=0`.
Likewise, nonzero shear makes `f''=2Pr*shear` incompatible with
`f(0)=f(infinity)=0`: a bounded solution with `f'(infinity)=0` would have
`f(infinity)=-2Pr*integral_0^infinity y*shear(y)dy < 0`. Thus even at
`Tw=1`, finite `Mr` does not produce the specified half-line thermal base.
The physical fixed-system-Mach path is `Mr=abs(Ro)*Msystem`; its exact
endpoint has `Mr=0` and requires an incompressible, not a `1/Ma^2`
compressible, stability operator.

## Solver behavior

`work/src/CompressibleBEKLinearStability.jl` now has
`NumericalParameters.thermal_model=:halfline` by default. It refuses
thermally active exact `Ro=0` cases, thermally active cases with no decaying
far-field heat mode (`Pr*Ro*W_inf>=0`), and cases whose source domain spans
fewer than ten thermal e-fold lengths. Explicit
`thermal_model=:finite_diagnostic` preserves those numerical experiments
but sets `evidence_scope` to `:ro0_thermal_degenerate`,
`:thermal_outflow_unverified`, or `:thermal_domain_unresolved`, as applicable.
The decay-length test is a numerical Gate, not a theorem about physical
existence. `PreparedCase` and output metadata record the far-field exponent
and number of decay lengths. `assemble_qep` rejects `Mr=0` with an explicit
incompressible-limit message. `at_fixed_system_mach` builds a matched path.

The public `work/scripts/run_compressible_bek_lsa_20260911.jl` CLI exposes
`--thermal_model=halfline|finite_diagnostic` and
`--min_thermal_decay_lengths=10`, and now uses the active operator-structure
and QEP directional-derivative Gates rather than the retired
compact/expanded audit. Its default output is under a new 2026-09-22 result
directory. Historical scripts from 2026-09-11/15 that silently ran
thermally singular finite-domain cases may require the explicit diagnostic
option to reproduce their old values; their recorded output remains
preserved.

## Numerical evidence

The final validation script
`work/scripts/validate_compressible_bek_thermal_limit_20260922.jl` completed
`8/8` guard tests and `24/24` algebraic/numerical tests. Selected cases:

| Ro | Mr | thermal length L | scope | alpha at R=116.3, beta=0.137 |
| ---: | ---: | ---: | --- | --- |
| -0.2 | 0.060 | 80 | derived candidate | `0.552239560-0.034492602i` |
| -0.1 | 0.030 | 160 | derived candidate | `0.540477199-0.018359507i` |
| -0.05 | 0.015 | 320 | derived candidate | `0.534940467-0.009511274i` |
| 0 | 0.100 | 40 | finite diagnostic only | `0.530157207-0.000748070i` |
| 1 | 0.100 | 40 | outflow diagnostic only | `0.461496796+0.141891420i` |

These are fixed-parameter spatial roots, not independently imposed neutral
conditions. The `Ro=-0.05,L=320` root changes by `4.30e-6` between
`N=49` and `N=59`; no full eigenvector-overlap continuation was done for
this comparison. For all tested QEPs, the named structural block residuals
are at most about `1.6e-15`, directional polynomial-derivative errors at
most `1.6e-9`, and relative eigenpair backward errors at most `2.7e-13`.
These tests check assembly/algebra and selected Ro/Co placements; they do
**not** independently prove every entry of the five-field operator.

The source-domain sweep in
`work/results/compressible_bek_thermal_limit_gate_20260922/low_ro_thermal_domain.csv`
shows why a simple sign test is insufficient. For `Ro=-0.05`, the fraction
`|f(0.75L)|/max|f|` falls from `0.1514` at `L=40` to `0.01127` at
`L=160` and `0.0001864` at `L=320`; `f(1)` changes from `-0.569592` to
`-0.576829`. The cases with approximately 11.5 thermal e-folds have
fractions below `2.23e-4` for all three sampled negative Ro values.
The ten-e-fold default is calibrated to this sweep and may need tightening
for other parameters.

An independent primitive-equation substitution was run on the admissible
`Ro=-1,L=40,N=120` and `Ro=-0.1,L=160,N=240` cases over `y<=12`. Maximum
absolute residuals for continuity are `1.66e-10/2.85e-10`, radial momentum
`2.82e-8/7.73e-9`, azimuthal momentum `9.65e-8/1.60e-7`, state `0/0`,
and energy `8.99e-9/3.61e-9`. Energy term-balance ratios are
`1.78e-6/2.61e-6`. These are local similarity boundary-layer residuals,
not full finite-radius Navier-Stokes residuals.

The repaired CLI was run for `Ro=-1,Mr=0.1,R=285.36,beta=0.07759,N=49`
and returned `alpha=0.385189621758-0.000392211915i`, backward error
`1.69e-14`, with all selected structure and boundary Gates near zero. This
is a von Karman endpoint smoke test, not a strict neutral-curve minimum.

## What is and is not validated

The thermal `Ro` placement, finite-domain labeling, far-field decay Gate,
selected operator structure, QEP algebra, low-Ro matched path and sampled
primitive base-flow residuals pass their stated tests. This does **not**
validate a finite-`Mr` compressible Ekman state at exact `Ro=0`, which the
assumed half-line similarity model cannot supply. It also does not resolve
the Bödewadt thermal outflow with a non-similar heat-transfer model, or
provide a fully independent five-equation linearized-operator Jacobian Gate
for every coefficient. Consequently, no claim is made that all `Ro`/`Mr`
neutral parameters are now physically correct.

Outputs: `work/results/compressible_bek_thermal_limit_gate_20260922/`;
source and runner changes are in `work/src/CompressibleBEKLinearStability.jl`
and `work/scripts/run_compressible_bek_lsa_20260911.jl`.

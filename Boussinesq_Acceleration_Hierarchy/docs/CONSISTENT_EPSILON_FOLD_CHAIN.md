# Epsilon-parameterised consistent-model fold chain

## Objective

Use

```text
epsilon = -Hinf > 0
```

as the continuation/control parameter for the acceleration-consistent BEK
model near the thermal-tail boundary.  The only requested production values
are

```text
epsilon = 0.1, 0.05, 0.02, 0.01, 0.005.
```

For every value, compute a directly corrected fold point
`(Ro_f(epsilon),Tw_f(epsilon))`, its profile and fixed-`Tw` Jacobian
diagnostics.  This study does not search for cusps, add further fixed-`Ro`
slices, derive asymptotics, or introduce stability and property variations.

## Augmented problem

For each prescribed `epsilon`, solve

```text
R(u,Tw;Ro) = 0,
J(u,Tw;Ro) v = 0,
v'v = 1,
Hinf(u) + epsilon = 0.
```

The unknowns are all collocation fields `u`, wall temperature `Tw`, the
right null vector `v`, and `Ro`.  The last equation makes `epsilon`, rather
than arclength or `Ro`, the local control parameter.  Solutions are followed
continuously from the verified `Ro=-0.5` fold; intermediate epsilon values
may be used internally by the corrector but are not added to the requested
result table.

## Numerical definition

- Fully self-consistent acceleration-density model in
  `work/src/BEKConsistent.jl`.
- `Pr=0.72`, `gamma=1`.
- Production map: `N=120`, `(a,b,c)=(2,0.6,0.5)`.
- Newton tolerance `1e-9`.
- The fixed-`Tw` Jacobian is row scaled before reporting its smallest
  singular value and singular-value ratio.
- Profiles are written only under a new directory in `work/results/`.
- Because the thermal length is `1/(Pr*epsilon)`, the two smallest epsilon
  values require a minimal `N=80,120,160` cross-check before any limiting
  claim.  Higher `N` is permitted only if that check has not converged.

## Outputs

- Shared Julia corrector: `work/src/BEKConsistent.jl`.
- Executable Julia study: `work/scripts/consistent_epsilon_fold_chain.jl`.
- Results: `work/results/consistent_epsilon_fold_chain/`.
- Requested clean table: `epsilon_fold_table.csv` in each run directory.

Conclusions and evidence scope must be appended to
`docs/RESEARCH_PROGRESS_LOG.md` after the computation is complete.

## Results

The direct fixed-`epsilon` fold corrector converged at all five requested
values.  The production `N=120` calculation is reliable through
`epsilon=0.02`; at the two smaller values it enters the resolution-sensitive
long-tail layer.  The recommended table therefore uses the common high-order
cluster `N=240,280,320`:

| `epsilon` | `Ro_f` | `Tw_f` | half-range in `Ro_f` | half-range in `Tw_f` |
|---:|---:|---:|---:|---:|
| 0.100 | -0.368218374 | 1.552811091 | `2.5e-12` | `5.0e-12` |
| 0.050 | -0.342987270 | 1.521896650 | `2.3e-10` | `2.8e-10` |
| 0.020 | -0.329321644 | 1.505245949 | `6.8e-8` | `8.3e-8` |
| 0.010 | -0.324984263 | 1.499970045 | `5.6e-6` | `6.8e-6` |
| 0.005 | -0.322827361 | 1.497347532 | `7.2e-5` | `8.7e-5` |

The increasing uncertainty at small `epsilon` is consistent with the thermal
length growing from about `13.9` at `epsilon=0.1` to `277.8` at
`epsilon=0.005`.  The highest-order temperature profiles remain monotone and
show this progressive tail extension directly.

Every point satisfies the augmented fold equations and has a row-scaled
fixed-`Tw` Jacobian singular value at roundoff (`O(1e-15)`).  Production
`N=120` residuals are below `4e-10`; the largest high-order corrector residual
is about `2e-9`.

## Finite-limit check

Both sequences become progressively flatter as `epsilon` is halved:

```text
epsilon:  0.020       0.010       0.005
Ro_f:    -0.329322   -0.324984   -0.322827
Tw_f:     1.505246    1.499970    1.497348
```

Thus the computed chain is consistent with finite limits for both
coordinates.  A linear intercept of only the last three computed points gives

```text
Ro_t diagnostic  = -0.320658670
Tw_t diagnostic  =  1.494709580.
```

This intercept is a finite-limit diagnostic, not a derived asymptotic law or
a final matched-asymptotic connection coordinate.  The clean numerical chain,
rather than the intercept, is the result to carry forward.

## Evidence files

- Recommended table:
  `work/results/consistent_epsilon_fold_chain/recommended_epsilon_fold_table.csv`.
- All resolution levels:
  `work/results/consistent_epsilon_fold_chain/epsilon_fold_convergence.csv`.
- Original production table:
  `work/results/consistent_epsilon_fold_chain/production_N120/epsilon_fold_table.csv`.
- Profiles and figures:
  `work/results/consistent_epsilon_fold_chain/epsilon_fold_profiles.png` and
  `recommended_epsilon_fold_chain.png`.
- Finite-limit diagnostic:
  `work/results/consistent_epsilon_fold_chain/finite_limit_diagnostic.txt`.

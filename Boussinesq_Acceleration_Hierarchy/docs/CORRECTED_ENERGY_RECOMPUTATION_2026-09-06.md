# Corrected Thermal Advection Recalculation

## Scope

The thermal equation was re-derived from the Lingwood BEK scaling. Since the
azimuthal velocity scale is `W*=ell*Omega*Ro*W` and the implementation stores
`W=-H`, the dimensionless energy equation is

```text
T'' + Pr*Ro*H*T' = 0.
```

The previous `T''-Pr*H*T'` form is identical only at `Ro=-1`; all non-von-
Karman thermal results using it are superseded. The corrected Jacobian,
`Ro` derivative, Hessian directional derivative and far-field length are
implemented in `work/src/BEKConsistent.jl`, `work/src/BEKThermal.jl` and
`work/src/BEKTraditionalForcing.jl`.

## Verification

The analytic derivatives were checked by finite differences:

| check | relative error |
|---|---:|
| consistent directional Hessian | `1.684876e-7` |
| consistent `Ro` residual derivative | `4.433743e-10` |

The von Karman traditional fold is unchanged within the spatial discretisation
error: `Tw=1.048021731`, `Hinf=-0.532753880` at degree 120 (the preserved
baseline is `Tw=1.048021731`, `Hinf≈-0.53276166`).

## Corrected topology at degree 80

The same continuation and classification criteria give:

| model | Ro=-1 | -0.75 | -0.5 | -0.25 | -0.1 |
|---|---|---|---|---|---|
| traditional | fold | thermal tail | thermal tail | thermal tail | thermal tail |
| consistent | smooth | smooth | fold | fold | thermal tail |

The consistent folds are at approximately `(Ro,Tw,Hinf)=(-0.5,1.69334,-0.08872)`
and `(-0.25,1.59222,-0.06205)`. Thus the qualitative fold-to-tail mechanism
survives, but its location is substantially different from the old,
uncorrected-energy calculation.

For `Ro>0` and the heated-wall branch with `Hinf<0`, the corrected far-field
thermal exponent has the wrong sign. A nontrivial localized semi-infinite
thermal solution is therefore not assigned in the corrected topology table;
this is an analytic localization obstruction, not a numerical failure.

## Corrected fold-tail chain

At `Ro=-0.5`, fixed-`Hinf` augmented folds converge for the larger tail scales.
Degree 280 gives:

```text
epsilon   Ro_f          Tw_f
0.10     -0.50427930   1.69816017
0.05     -0.48619670   1.67818717
0.02     -0.47642533   1.66784029
0.01     -0.47331991   1.66461919
```

The `epsilon=0.005` point remains resolution-sensitive and is not used for an
independent asymptotic coefficient. The corrected outer slow-manifold solve
is converged by `Ymax=60` and gives

```text
Ro_t = -0.4703674163,   Tw_t = 1.6615874195.
```

This supports a finite fold-to-tail endpoint and the same thermo-rotational
slow-tail interpretation, but the old limiting coordinates and old linear
coefficients must be discarded. No claim of a new coefficient is made here.

## Evidence boundary

Outputs are under `work/results/corrected_energy_*`; preserved baseline files
were not overwritten. This report does not modify the learning manual.

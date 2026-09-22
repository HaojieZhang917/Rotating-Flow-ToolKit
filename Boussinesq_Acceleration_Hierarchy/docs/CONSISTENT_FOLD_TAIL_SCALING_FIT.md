# Scaling fit for the consistent fold--tail limit

## Scope

Use only the existing recommended fixed-`epsilon` fold data at

```text
epsilon = 0.1, 0.05, 0.02, 0.01, 0.005
```

to test

```text
Ro_f(epsilon) - Ro_t ~ A_R epsilon^p,
Tw_f(epsilon) - Tw_t ~ A_T epsilon^q.
```

No new base-flow solve, `Ro` slice, smaller `epsilon`, stability calculation
or matched-asymptotic derivation is part of this task.

## Input dependency

The sole numerical input is

`work/results/consistent_epsilon_fold_chain/recommended_epsilon_fold_table.csv`,

which records the `N=240--320` consensus values and their high-order
half-ranges.  Its model definition and convergence scope are documented in
`docs/CONSISTENT_EPSILON_FOLD_CHAIN.md` and log entry L011.

## Tests

1. Compute the limit-independent local exponent from the final geometric
   triple `epsilon=(0.02,0.01,0.005)`:

   ```text
   p_eff = log((y(0.02)-y(0.01))/(y(0.01)-y(0.005)))/log(2).
   ```

2. Fit both

   ```text
   y = y_t + A epsilon
   y = y_t + A epsilon + B epsilon^2
   ```

   on all five points and after removing the largest point `epsilon=0.1`.
   Compare the intercept and linear coefficient across windows.

3. For the selected small-epsilon quadratic fit, check that
   `(y-y_t)/epsilon` tends to `A` and that
   `(y-y_t-A epsilon)/epsilon^2` is approximately constant.

4. Report residual norms without interpreting the fitted polynomial as a
   derived asymptotic expansion.  The purpose is to identify the distinguished
   scaling required by a later analysis.

## Outputs

- Julia analysis: `work/scripts/fit_consistent_fold_tail_scaling.jl`.
- Results: `work/results/consistent_fold_tail_scaling/`.
- Conclusions must be appended to `docs/RESEARCH_PROGRESS_LOG.md` before use
  in a manuscript claim.

## Results

### Exponents

The final geometric triple `epsilon=(0.02,0.01,0.005)` gives the
limit-independent difference ratios

```text
Delta Ro(0.02->0.01) / Delta Ro(0.01->0.005) = 2.01093,
Delta Tw(0.02->0.01) / Delta Tw(0.01->0.005) = 2.01177.
```

Hence

```text
p_eff = 1.007864,
q_eff = 1.008468.
```

Perturbing every input independently over its documented high-order
half-range gives

```text
p_eff in [0.9551,1.0625],
q_eff in [0.9557,1.0631].
```

An independent unknown-limit power fit over the four smallest points gives
`p=1.0467`, `q=1.0484`.  Both diagnostics are consistent with, and do not
resolve a meaningful departure from, unit exponent.  The supported numerical
scaling is therefore `p=q=1`.

### Polynomial fits and window stability

| variable | fit window | limit | linear coefficient | quadratic coefficient |
|---|---|---:|---:|---:|
| `Ro_f` | all five | -0.320760492 | -0.415085998 | -0.594667305 |
| `Ro_f` | omit `epsilon=0.1` | -0.320719730 | -0.420342495 | -0.500001415 |
| `Tw_f` | all five | 1.494843888 | 0.503275261 | 0.763601163 |
| `Tw_f` | omit `epsilon=0.1` | 1.494786622 | 0.510660114 | 0.630605042 |

Removing the largest `epsilon` changes the intercepts by only
`4.1e-5` (`Ro`) and `5.7e-5` (`Tw`), and changes the linear coefficients by
about `1.3--1.5%`.  Direct corner perturbation of the four-point input gives

```text
A_R in [-0.42766,-0.41303],
A_T in [ 0.50178, 0.51955].
```

A concise coefficient estimate is therefore

```text
A_R approximately -0.42,
A_T approximately  0.51.
```

The preferred four-smallest-point representation is

```text
Ro_f = -0.32071973 - 0.42034250 epsilon
       - 0.50000142 epsilon^2 + residual,

Tw_f =  1.49478662 + 0.51066011 epsilon
       + 0.63060504 epsilon^2 + residual.
```

### Residual-order check

For all five points, adding the quadratic term reduces RMS residuals from
`5.09e-4` to `1.68e-5` for `Ro` and from `6.54e-4` to `2.32e-5` for `Tw`,
approximately factors of `30` and `28`.  After omitting `epsilon=0.1`, the
reductions are factors of about `13` for both variables.

With the four-point quadratic fit, `(y-y_t)/epsilon` tends toward the fitted
linear coefficient.  Over the resolved intermediate points,
`(y-y_t-A epsilon)/epsilon^2` is close to the fitted `B`; its scatter at
`epsilon=0.01,0.005` is consistent with the documented amplification of input
uncertainty after division by `epsilon^2`.

Thus the numerical residual structure supports an `O(epsilon)` leading
correction and an approximately `O(epsilon^2)` next correction.  It does not
derive the coefficients analytically.

## Supported distinguished scaling

The result to pass to a later matched-asymptotic analysis is

```text
Ro - Ro_t = O(-Hinf),
Tw - Tw_t = O(-Hinf),
```

with numerical first-order coefficients approximately `-0.42` and `+0.51`.

Evidence files are under
`work/results/consistent_fold_tail_scaling/`, particularly
`polynomial_fits.csv`, `exponent_diagnostics.csv`,
`uncertainty_sensitivity.csv`, `normalized_scaling.csv`,
`scaling_fits.png` and `scaling_residuals.png`.

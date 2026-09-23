# Legacy `sim`-coordinate `Ro=0` grid and far-field diagnostic (2026-09-22)

## Scope and provenance

This is a numerical diagnostic of the *current working-tree* legacy notebook
chain, not a validation of the corrected general-`Ro` BEK equations. We keep
`Ro=0`, `Co=2`, `Tw=1`, `Mr=0.1`, `R=116.3`, `beta=0.137`, `omega=0`,
`gamma=1.4`, `Pr=0.72`, `Ma=Mr/R`, the old `baseflow_var`, `T_ca`,
`Spatial_mode_BEK`, `assemble_mat`, boundary deletion, PEP, and IAR. Profile
interpolation uses `interp(...,"sim")`. A single legacy base-flow BVP solution
is reused for all stability runs. The source files were already modified in
the user's working tree before this study and were not edited here:

- `compress/BEK/roughness/CRD_STA.jl` SHA-256
  `57d8f2cc4a94ba9a6c7431aebd12f37d308ff48f45da451ebf2a8b2fadc9084e`
- `compress/BEK/roughness/CRD_STA_CACU.ipynb` SHA-256
  `5e4f3d0603cfd6104bd64032605c4e2b7913c71773cae2350f572a2e59d9df45`

The executable is `work/scripts/check_legacy_simcoord_ro0_grid_farfield_20260922.jl`.
It uses `iar` with target `0.5`, four returned modes, `maxit=500`, and
`tol=1e-12`. We identify the reference root as the one nearest the notebook's
saved `0.5403618717476425+0.00033737176837288875i`; all other cases select
the mode with greatest normalized five-field eigenvector overlap on the common
near-wall interval `0<=y<=8`. Reported backward errors are the relative QEP
residuals. The baseline custom grid reproduces legacy `Cheb(69)` coordinates
exactly and `D`, `D2` to floating-point roundoff (maximum differences
`2.84e-14`, `1.46e-11`).

## Results

| Variation | Parameters | Tracked alpha | Overlap |
| --- | --- | --- | ---: |
| Reference | `N=69,a=6,cap=20` | `0.540361211263+0.000337288886i` | 1 |
| Grid | `N=49` | `0.540358639479+0.000344005170i` | 0.999799 |
| Grid | `N=59` | `0.540361658197+0.000337521706i` | 0.999917 |
| Grid | `N=79` | `0.540361386082+0.000337126261i` | 0.999953 |
| Grid | `N=89` | `0.540361335286+0.000337472556i` | 0.999953 |
| Map scale | `N=69,a=4,cap=20` | `0.540361240582+0.000337493185i` | 0.999956 |
| Map scale | `N=69,a=8,cap=20` | `0.540361347761+0.000337468110i` | 0.999935 |
| Map scale | `N=69,a=10,cap=20` | `0.540361210363+0.000336705935i` | 0.999888 |
| Sampling cap | `N=69,a=6,cap=30` | `0.540361656843+0.000337300691i` | ~1 |
| Sampling cap | `N=69,a=6,cap=40` | `0.540361446653+0.000337361062i` | ~1 |

The full-precision data, relative QEP backward errors (`9.98e-16` to
`1.52e-13`), and counts of duplicated capped coordinate nodes are in
`work/results/legacy_simcoord_ro0_convergence_20260922/modes.csv`.

## Conclusion and limitations

The *selected legacy eigenbranch* is numerically stable under these grid and
far-field sampling variations. Its imaginary part remains about
`3.37e-4`, so it satisfies the user-chosen approximate-neutral criterion
`|Im(alpha)|<1e-3`; this is not an exact-neutral or critical-curve calculation.
The reproduced reference differs from the notebook's saved value by
`6.66e-7`, consistent with a small solver/environment variation; the run is
tied to the source hashes above.

The legacy rational Chebyshev map has an infinite endpoint but clips sampled
coordinates at `cap` while leaving its derivative matrices unchanged. At
`N=69,a=6,cap=20`, 22 nodes sample the same `y=20`; at caps 30 and 40 there
are 18 and 15 such nodes, respectively. Changing the cap is therefore a
**sampling-sensitivity diagnostic**, not a mathematically consistent finite
domain-length convergence test. Varying `a` tests the rational-map resolution
with this same legacy clipping. These numerical checks do not establish the
physical consistency of the old heat BVP/operator, the singular `Ro=0`
thermal limit, the `Mr=0.3` model, other `Ro` values, or the minimum of a
neutral curve in `(R,beta)`.

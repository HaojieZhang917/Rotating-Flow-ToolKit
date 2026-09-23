# `Ro=0` legacy/project comparison at matched notebook parameters

Date: 2026-09-22

## Question and controlled setting

The reported inability of the project solver to find an Ekman near-neutral
mode was checked at the *notebook's own parameters*, not at a different Mach
number or spanwise wavenumber: `Ro=0`, `Co=2`, `Mr=0.1`, `Tw=1`, `R=116.3`,
`beta=0.137`, `omega=0`, `gamma=1.4`, `Pr=0.72`, `N=69`. The project solver uses
its rational half-line grid, its corrected general-`Ro` candidate operator,
and its current thermal BVP. The targeted QEP used `target=0.54`, `neigs=6`,
`tol=1e-11`, balanced matrices. The executable is
`work/scripts/compare_ro0_matched_notebook_parameters_20260922.jl`; the full
mode list is in `work/results/ro0_matched_notebook_comparison_20260922/modes.csv`.
The project frontend SHA-256 was
`4ec9159275f50162b80df4529cbbe8208997f7186cb66f21840e887f00fa072c`;
the preserved backend SHA-256 was
`399e56eb6f9581d6da3e9d19749d7146d64dcffacae252a9ba7c328dab7eb76d`.
As in earlier studies, the script replaces the broken
`using DifferentialEquations` meta-package import in memory with
`using BoundaryValueDiffEq`; it does not edit either solver source.

## Result

| Solver | Matched-point alpha | `|Im(alpha)|` | Relative QEP backward error |
| --- | ---: | ---: | ---: |
| Legacy `sim` notebook chain | `0.540361211263+0.000337288886i` | `3.37e-4` | `1.51e-13` |
| Project corrected-BEK candidate | `0.530157206567-0.000748070374i` | `7.48e-4` | `7.45e-15` |

Both meet the user's `|Im(alpha)|<1e-3` *approximate-neutral* criterion at
these parameters. The two complex roots differ by `0.010261565` and neither
result alone establishes that their eigenvectors represent an identical
physical branch across the two discretized operators. The project temperature
varies by `max|T-1|=0.00134293228`, and its case is explicitly labeled
`ro0_thermal_degenerate`.

## Why the roots need not coincide

This is not a same-equation comparison. At `Ro=0`, the legacy `f_q` still solves
`f''=2Pr[(u_old')^2+(v_old')^2+u_old f-phi_old f']`, with the pre-flip legacy
source profiles, whereas the project's unified equation
reduces to `f''=2Pr[(U')^2+(V')^2`. Even with `Tw=1`, `T=1-(gamma-1)Mr^2 f/2`,
so a change in `f` changes the compressible base state. The old stability
operator also leaves its cylindrical continuity curvature (`D_11=rho`) and
base normal transport (`C_14=rho H`) active at `Ro=0`; the corrected candidate
multiplies those terms, and the associated energy normal-transport terms, by
`Ro`, so they vanish in the Ekman limit. The grids and far-field profile
sampling also differ; the legacy grid duplicates clipped far-field sample
positions, whereas the project grid retains a unique `Inf` endpoint.

The current check identifies **structural differences**, not a measured
contribution of each difference to the `0.01026` root gap. Such attribution
would require a one-change-at-a-time operator/base-state cross-evaluation and
cross-model eigenvector matching. Prior project calculations at `Mr=0.3`,
`beta=0.132` are at different parameters and cannot be used to claim that
the matched `Mr=0.1`, `beta=0.137` point has no near-neutral root.

## Evidence boundary

The project case's `ro0_thermal_degenerate` label remains important: the
unified `Ro=0` heat equation loses normal transport, and a finite-domain
thermal solution with nonzero `Mr` is not yet a validated half-line base
state. A near-neutral QEP root, in either model, does not settle that physical
issue or establish the minimum of a neutral curve. No claim is made here
about `Mr=0.3`, `Ro=1`, or full three-dimensional stability.

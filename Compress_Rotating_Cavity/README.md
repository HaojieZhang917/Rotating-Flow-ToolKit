# Compressible rotating-cavity project scaffold

This directory now contains a reproducible isothermal benchmark pipeline. It
is the baseline that will be extended to the thermal/Boussinesq and
compressible models.

## Program layout

- `BaseFlow_cavity.jl`: original notebook-compatible base-flow entry point.
- `receptivity.ipynb`: original exploratory notebook and reported reference values.
- `src/CRC_BF.jl`: organized copy of the legacy SciPy/PyCall base-flow solver.
- `src/CRC_STA.jl`: the four-field incompressible/isothermal local spatial
  stability operator used in the paper and retained as the regression
  reference.
- `src/TraditionalBoussinesq.jl`: traditional centrifugal Boussinesq thermal
  similarity equations and pressure continuation.
- `src/BlackburnBoussinesq.jl`: acceleration-consistent Blackburn closure,
  including the `chi=2-T` mass and momentum weighting.
- `src/TraditionalBoussinesqTemporal.jl`: common axisymmetric,
  similarity-preserving temporal stability operator. Use
  `closure=:traditional` or `closure=:blackburn`.
- `src/TraditionalBoussinesqLST.jl`: five-field local perturbation operator
  `(u,v,w,theta,p)` for the traditional Boussinesq closure. Its isothermal
  `[u,v,w,p]` block is checked against `CRC_STA.jl` before boundary reduction
  and after boundary/gauge elimination.
- `src/CompressRotatingCavity.jl`: Julia API joining the two solvers.
- `scripts/benchmark_isothermal.jl`: command-line regression benchmark.
- `scripts/benchmark_traditional_boussinesq.jl`: complete Boussinesq fold
  benchmark.
- `scripts/benchmark_temporal_stability.jl`: three-branch temporal spectrum
  benchmark at `Tw=1.16`.
- `scripts/benchmark_traditional_lst.jl`: traditional thermal LST module and
  isothermal regression against the paper operator.
- `scripts/verify_paper_isothermal.jl`: checks the manuscript's isothermal
  spatial benchmark, the `Re_h=1000` core parameters, and the rounded temporal
  critical points in Table IV.
- `scripts/compare_closures.jl`: same-grid comparison of the traditional and
  Blackburn closures, including branch status, leading growth rate, and
  normalized mode-shape tables.
- `scripts/scan_traditional_lst_temperature.jl`: fixed-parameter local temporal
  LST scan along the lower traditional branch, using the paper's Type-I and
  Type-II isothermal reference points. The same run also tracks the spatial
  PEP root `alpha` at fixed `omega` and `beta`.
- `scripts/verify_lower_branch_convergence.jl`: continuation and grid
  convergence of the isothermal-connected lower/principal branch.
- `data/reference/`: archived numerical reference profile used only to start
  the isothermal branch.
- `benchmark_results/`: summary and branch data produced by the benchmarks.

The legacy `CRC_BF`/`CRC_STA` path remains isothermal and incompressible. In
that path, mode 1 uses `Ts` as the axial wall mass-flux parameter in
`H(0)=-Ts`; it is not a temperature. The thermal modules listed above use
`Tw` for the rotor wall temperature.

## Run

```bash
julia --project=. scripts/benchmark_isothermal.jl
```

The benchmark targets the values recorded in cell 7 of `receptivity.ipynb`.
It checks the wall boundary conditions, direct eigenvalue, transpose-adjoint
eigenvalue, and receptivity coefficient.

To verify the isothermal values reported in the manuscript, run:

```bash
julia --project=. scripts/verify_paper_isothermal.jl
```

The Table-I/Fig.-16 check uses the spatial polynomial problem. Table-IV values
are checked in the temporal generalized eigenproblem at the reported real
`(R, alpha, beta)` values; the table is rounded, so the computed imaginary
frequency is reported as a neutrality residual.

## Traditional Boussinesq benchmark

The second benchmark retains the radial pressure gradient as an unknown and
continues in that parameter. This is required to pass through the saddle-node
where `dTw/dPi=0`.

```bash
julia --project=. scripts/benchmark_traditional_boussinesq.jl
```

For `Re_h=1000` and `Pr=0.72`, the principal fold is reproduced at
`Tw=1.167676484...`. The branch and summary are written to
`benchmark_results/`.

## Temporal stability and Blackburn comparison

```bash
julia --project=. scripts/benchmark_temporal_stability.jl
julia --project=. scripts/compare_closures.jl
julia --project=. scripts/scan_traditional_lst_temperature.jl
```

The temporal calculation is restricted to the axisymmetric similarity-
preserving subspace used by the saddle-node dynamical report; it is not the
full finite-radius three-dimensional spectrum. At `Tw=1.16`, the traditional
closure gives leading growth rates `(-0.0058285, +0.0028491, -0.0042999)` for
the upper, middle, and principal branches. In the Blackburn comparison, the
principal branch has `-0.0277213` and the second converged branch has
`-0.0154466 ± 0.1438453 i`. The traditional principal fold is retained at
`Tw=1.167676484...`; no Blackburn principal fold was detected in the
continuation scan up to `Tw=1.75106` (the linear factor reaches its formal
wall degeneracy at `chi_w=2-Tw=0` when `Tw=2`).

The comparison CSV files also record the normalized radial, azimuthal,
axial, and temperature mode components. A failed Newton seed is retained in
the status table rather than being silently reported as a physical branch.
`closure_critical_comparison_Re1000.csv` evaluates both closures at the
traditional fold temperature: the traditional leading eigenvalue is
approximately zero, while the Blackburn state remains damped.

The temperature-scan CSV filters numerical near-infinite eigenvalues generated
by the pressure-null mass matrix before ordering the temporal spectrum. It
contains both the most amplified finite mode and the continuously tracked mode
near each manuscript reference frequency.

## Lower/principal branch verification

```bash
julia --project=. scripts/verify_lower_branch_convergence.jl
```

This is the primary stability workflow. It follows the branch connected to
the isothermal state at `Tw=1.16, 1.40, 1.60, 1.75` and repeats the
calculation at Chebyshev degrees 60, 80, and 100. The Blackburn branch remains
temporally damped over this scan; its leading eigenvalue becomes a damped
oscillatory pair above the low-temperature range. The traditional fold
temperature is converged to the displayed digits across the three grids.

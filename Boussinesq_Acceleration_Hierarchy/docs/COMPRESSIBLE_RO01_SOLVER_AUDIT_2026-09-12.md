# Compressible `Ro=0,1` solver audit (2026-09-12)

## Registered scope before calculation

This diagnostic follows the legacy similarity-coordinate comparison at
`Mr=0.3`, `Tw=1`, `N=69`.  It is intended to separate four questions that
cannot be answered by comparing the first eigenvalue returned by IAR:

1. whether the corrected thermal base state has the expected `Ro=0` limit;
2. whether the inherited `Ro=-1 -> +1` operator branch is continuous;
3. which spatial mode is obtained by eigenvector-overlap continuation from a
   finite low-Mach seed to `Mr=0.3`;
4. whether the tracked endpoint is neutral at the supplied fixed `(R,beta)`.

The calculation is a frozen-local-base-flow, spatial linear-stability
diagnostic.  It is not a self-consistent neutral-curve continuation, a full
three-dimensional stability result, or a validation of the arbitrary-`Ro`
compressible equations.

## Model and parameters

The diagnostic uses `work/src/CompressibleBEKLinearStability.jl` and the
operator implementation in
`/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/src/CRD_STA.jl`.
The base state uses the corrected two-point thermal BVP and the unified
similarity-coordinate convention in the former file.  The stability operator
uses the explicitly labelled `unverified_bek_extension` for `Ro=0,1`.

Fixed parameters are

- `Tw=1`, `gamma=1.4`, `Pr=0.72`, `omega=0`;
- Ekman: `Ro=0`, `R=116.25`, `beta=0.132`, supplied `alpha_ref=0.2775`;
- Bodewadt: `Ro=1`, `R=27.26`, `beta=0.127`, supplied
  `alpha_ref=0.5231`;
- semi-infinite rational similarity grid, with the endpoint represented by
  `y=Inf`, finite base-flow BVP support through `y=40`, and production check
  at `N=69`;
- Mach path `Mr=(0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30)`.

At the first Mach point the targeted eigensolve is seeded by `alpha_ref`.
Every later solve is seeded by the preceding eigenvalue and eigenvector, and
the selected mode is ranked first by normalized eigenvector overlap.  The
polynomial backward error and overlap are retained at every step.

## Source-data and runtime dependencies

- supplied reference tuple and the legacy comparison:
  `docs/LEGACY_SIMCOORD_RO01_COMPARISON_2026-09-11.md`;
- inherited BEK base-flow and compressible QEP source:
  `/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/src/CRD_STA.jl`;
- Turkyilmazoglu & Uygun (2006), equations (5)--(15), only for the strict
  rotating-disk (`Ro=-1`) slice;
- Lingwood (1997) for the incompressible BEK family and its `Ro` convention.

The current root Julia environment cannot load the `DifferentialEquations`
meta-package because its installed `OrdinaryDiffEqBDF` expects the missing
`SciMLBase.AutoDespecialize` symbol.  The audit script therefore evaluates an
in-memory copy of the same `CRD_STA.jl` source with only the unused
`using DifferentialEquations` imports replaced by
`using BoundaryValueDiffEq`.  No original source or package manifest is
modified.  This workaround is acceptable only for this diagnostic and is
recorded in the output metadata.

The inherited `Cheb` routine formed its differentiation matrices using an
algebraic map to infinity but clipped every reported coordinate above 40 to
40.  This mixed representation produced repeated far-field coordinates and
made cross-grid eigenvector interpolation ambiguous.  The registered solver
now retains the same algebraic differentiation map, reports the terminal node
as `Inf`, and assigns the exact asymptotic base state to nodes beyond the
finite BVP support.  Grid-to-grid mode transfer is performed in the compact
coordinate `s=y/(1+y)`.

## Results

### Semi-infinite grid representation

The rational grid now has one explicit endpoint at `y=Inf`; all preceding
coordinates are finite and strictly increasing.  For `N=59,69,79,89,99`, the
old clipped representation contained respectively 13, 15, 17, 20 and 22
copies of `y=40`.  The corresponding new grids contain exactly `N` finite
nodes plus the one infinite endpoint.  The differentiation matrices are
identical to the inherited matrices to machine precision (reported relative
changes in both `D` and `D2` are exactly zero), because only the inconsistent
coordinate reporting and base-state sampling were corrected.

Nodes beyond the base-flow BVP support at `y=40` are assigned
`F=0`, `G=1`, `T=rho=1` and the computed constant `H_inf`.  Cross-grid mode
transfer uses the compact coordinate `s=y/(1+y)`, so the infinite endpoint is
represented exactly by `s=1` without extrapolation.

### Corrected thermal base state

At `Mr=0.3`, the `Ro=0` energy equation gives an isothermal base to roundoff:
`max|T-1|=3.33e-16`.  At `Ro=1`, the corrected solution has
`max(T)=1.1226852428` and `max|rho-1|=0.1092783962`.  Thermal BVP endpoint
errors are below `3.5e-17`.

### Mach continuation at `N=69`

Eigenvector-overlap continuation along
`Mr=(0.01,0.03,0.05,0.10,0.15,0.20,0.25,0.30)` gives at the final point

- Ekman (`Ro=0`):
  `alpha=0.315366030307-0.088067397220i`, final overlap effectively 1;
- Bodewadt (`Ro=1`):
  `alpha=0.563467793123-0.150518845867i`, final overlap 0.998749.

Both branches remain plainly non-neutral at the supplied fixed `(R,beta)`.
Thus first-root ordering is not the main explanation for the discrepancy on
these reference-targeted branches.

### Grid check at `Mr=0.3`

For Ekman, `N=59,69,79` give respectively

```text
0.315366030136-0.088067397556i
0.315366030307-0.088067397220i
0.315366030304-0.088067397218i
```

so this candidate is grid converged.  For Bodewadt, direct target solves at
`N=59,69,79` give

```text
0.552452550372-0.179782289433i
0.563467756147-0.150518846288i
0.567409336863-0.142941057878i
```

and remain grid sensitive.  Cross-grid eigenvector continuation through
`N=59,69,79,89,99` has overlaps
`seed, 0.993996, 0.999675, 0.997254, 0.999085`, but the endpoint sequence still
oscillates and reaches
`0.559223368023-0.160621751355i` at `N=99`.  The high overlaps show that this
is not merely an arbitrary nearest-root switch; the rational-grid sequence is
not yet converged for the Bodewadt candidate.

### Conclusion and evidence scope

Mapping the grid consistently to infinity removes a real coordinate/sampling
defect, but it leaves the tracked eigenvalues essentially unchanged and does
not resolve the `Ro=0,1` mismatch.  The remaining investigation must still
address the unverified general-`Ro` operator/sign conventions and, for
Bodewadt, the unresolved grid-sensitive eigenvalue cluster.  These are
frozen-local-base-flow, similarity-subspace spatial-stability diagnostics,
not self-consistent neutral curves or full three-dimensional stability
results.

Numerical evidence is stored in
`work/results/compressible_solver_ro01_audit_20260912/`, especially
`semi_infinite_grid_mapping.csv`, `mach_mode_tracking_N69.csv`,
`endpoint_grid_convergence_Mr03.csv`, and
`bodewadt_grid_overlap_continuation_Mr03.csv`.

# Corrected compressible BEK fixed-beta neutral study (2026-09-12)

## Registered calculation scope

This study tests the corrected general-`Ro` compressible BEK operator at the
two supplied stationary-wave reference tuples.  It seeks a zero of
`imag(alpha)` by continuous spatial-mode tracking at fixed `beta`; it does not
claim to locate the nose of a two-parameter neutral curve.

The calculations are frozen-local-base-flow, similarity-subspace spatial
stability calculations.  They are not full three-dimensional stability
calculations.  The `Ro=0,1` equations are a derived BEK candidate whose strict
paper validation slice is `Ro=-1`.

## Model and fixed parameters

- solver: `work/src/CompressibleBEKLinearStability.jl`;
- legacy backend, read without modification:
  `/home/zhj/Rotating-Flow-ToolKit/RotatingDiskFlow/src/CRD_STA.jl`;
- similarity coordinate: semi-infinite rational Chebyshev map, with the last
  node represented explicitly by `Inf`;
- `Tw=1`, `Mr=0.3`, `omega=0`, `gamma=1.4`, `Pr=0.72`;
- Ekman: `Ro=0`, fixed `beta=0.132`, supplied reference `R=116.25`,
  `alpha_ref=0.2775`;
- Bödewadt: `Ro=1`, fixed `beta=0.127`, supplied reference `R=27.26`,
  `alpha_ref=0.5231`.

The reference-point mode is first followed on the Mach path
`Mr=(0.01,0.05,0.10,0.15,0.20,0.25,0.30)`.  Its eigenvector and eigenvalue are
then used to continue in both directions in `R` at `N=49`.  The registered
coarse ranges are `R=40..250` for Ekman and `R=10..75` for Bödewadt.  Any
adjacent sign change of `imag(alpha)` is refined with a safeguarded secant
iteration.  A located zero is re-solved at `N=59` and `N=69` using the
continued eigenvalue as the target.

At every point the selected eigenvalue, normalized eigenvector overlap with
the preceding point, step size, and polynomial backward error are retained.
Failure to find a sign change means only that the tracked branch has no fixed-
`beta` crossing over the registered range.

## Equation/source dependencies

- Turkyilmazoglu & Uygun (2006), equations (5)--(15), validates only the
  compressible von Karman (`Ro=-1`) endpoint;
- William's thesis, chapter 2 for the compressible rotating similarity base
  state and chapter 3 as a perturbation-equation consistency check;
- Lingwood (1997) general-BEK convention `Co=2-Ro-Ro^2` and the local
  incompressible placement of `Ro` and `Co`;
- operator audit: `docs/COMPRESSIBLE_BEK_RO_CO_OPERATOR_AUDIT_2026-09-12.md`.

The root Julia environment currently cannot load the
`DifferentialEquations` meta-package because of an installed SciML version
mismatch.  The executable study therefore evaluates in memory the same source
after replacing only the unused import with `BoundaryValueDiffEq`; neither
the preserved backend nor a package manifest is edited.

Before interpreting the supplied tuples, a separate low-Mach endpoint gate is
also evaluated at the rounded Lingwood (1997) stationary-wave nose values:
`(Ro,R,beta,alpha)=(-1,290.1,0.077,0.381)`,
`(0,116.3,0.137,0.528)`, and `(1,27.4,0.115,0.487)`, at `Mr=0.01` and
`N=(49,59)`.  These rounded incompressible values are a physics gate, not exact
compressible neutral roots.

A second, reference-anchored continuation is registered because targeting the
supplied `alpha_ref=(0.2775,0.5231)` can select a different spatial branch.
For each flow it starts from the corresponding low-Mach Lingwood mode, follows
`Mr` to `0.3`, then follows `beta` to the supplied value, and only then varies
`R` to impose `imag(alpha)=0`.  This result is the physically anchored fixed-
`beta` neutral candidate; the supplied-alpha scan is retained as a separate
mode diagnostic.

For grid convergence, neutrality is re-imposed independently at each of
`N=(49,59,69)` inside the brackets `R=(116.3,125)` for Ekman and `R=(30,35)`
for Bödewadt.  Merely evaluating a finer grid at the coarse-grid neutral `R`
is retained only as a residual check, not reported as the refined neutral
parameter.

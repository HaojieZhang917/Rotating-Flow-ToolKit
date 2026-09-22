# Research execution plan

## Priority override: finite-radius physical validation (2026-09-09)

The immediate priority is no longer further refinement of the
similarity-space matched-adjoint fold direction. The primary question is
whether the similarity fold--tail topology survives in a radially
non-similar heated BEK flow.

The required order is now:

1. derive and regression-test the non-similar axisymmetric equations;
2. compute the linear radial transverse pencil and forced response of the
   similarity subspace;
3. specify a physical symmetry breaking and its radial/edge boundary data;
4. determine whether a parabolized radial formulation is well posed;
5. solve three benchmark Rossby cases in a global finite-radius model;
6. resume matched-adjoint work only if the global solution develops a
   similarity-like interior region near the fold--tail limit.

The first two diagnostics are recorded in
'docs/NONSIMILAR_HEATED_BEK_BOUNDARY_LAYER_2026-09-09.md'. They establish
that similarity is an exact invariant subspace and that the existing
generic-tail radial velocity changes sign with wall-normal position for
\(Ro=-0.50\) and \(-0.40\), while its first-order part degenerates near
\(Ro_t\). Consequently, a one-way radial march is not the primary physical
solver. Before that global calculation, the pencil
\((J_{\rm sim}+\lambda M_r)\widehat Q=0\) will quantify spatial radial
sensitivity of the similarity invariant subspace. The preferred final
validation is a global axisymmetric finite-disk
calculation including radial diffusion and the edge region, ideally using
the existing compressible solver under matched rotation, thermal, property,
and geometry definitions.

Until that calculation is complete, the similarity fold--tail and
matched-adjoint results must be described as structures of the
self-similar invariant subspace, with their finite-radius physical
persistence unresolved.

## Research question

Determine how the BEK differential-rotation Rossby parameter changes the error introduced by Boussinesq density--acceleration closures, and how that error is amplified into changes in base-flow saddle-nodes, branch topology and stability predictions.

The intended causal chain is

```text
BEK differential rotation
    -> acceleration and base-flow structure
    -> traditional versus canonical closure residual
    -> fold/cusp displacement and critical amplification
    -> similarity-subspace and Type-I/Type-II stability errors.
```

`Ro` is the physical BEK family parameter. Term coefficients `eta_j` are diagnostic sensitivity/homotopy coordinates, not physical parameters and not empirical coefficients to be optimised.

## Gate 1: dimensional derivation and BEK reduction

Derive the dimensional rotating-frame momentum equation and separate local, convective/streamline-curvature, Coriolis, frame-centrifugal, Euler and imposed-body accelerations. Reproduce Lingwood's BEK definitions and obtain the similarity equations with

```text
Ro=-1: von Karman,
Ro=0:  Ekman,
Ro=1:  Bodewadt,
Co=2-Ro-Ro^2.
```

Derive the traditional and canonical thermal equations from the same parent formulation. Identify which gradient accelerations may be absorbed into the reference pressure before multiplying by density fluctuation, and verify reference-frame consistency.

Deliverables:

- checked notation, dimensional table and sign convention;
- derivation of the steady BEK ODEs and boundary conditions;
- term-by-term map from dimensional acceleration to the reduced equations;
- exact endpoint-model definitions and the permitted diagnostic `eta_j` paths;
- derivation and reference-frame checks recorded in `docs/` before numerical claims.

Stop condition: do not implement a parameter sweep until the three BEK limits and the traditional/canonical endpoint map are unambiguous.

## Gate 2: verified self-consistent steady models

Implement reusable Julia code in `work/src/` for one parent steady residual with physically derived endpoint configurations. Use `eta_j` only for derivative checks and diagnostic continuation.

Verification sequence:

1. recover the isothermal BEK profiles across representative `Ro` values;
2. recover the von Karman, Ekman and Bodewadt limits;
3. reproduce the preserved traditional von Karman fold at `Ro=-1`;
4. confirm boundary residual, spectral-order and far-field convergence;
5. compare endpoint profiles at `Tw=1` and selected small temperature differences.

Deliverables:

- Julia module and unit/regression tests;
- an isothermal BEK benchmark table;
- endpoint base-flow profiles and residual audits;
- a documented production-sweep configuration under `docs/`;
- all new outputs under a new directory in `work/results/`.

## Priority Gate R3: radial transversality of the fold--tail limit

Before starting a global finite-disk calculation, test the steady
axisymmetric radial-response pencil along the corrected finite-ϵ fold
chain. Use the normalization-invariant coefficient

```text
chat_c_r = <w,M_r v>/<w,R_Tw>
```

rather than a raw coefficient based on a potentially singular `w'v=1`
normalization. Validate the fold normal form with direct near-fold pencil
eigenpairs, and determine the power of ϵ with grid and mapping checks.

Decision: if `chat_c_r` tends to a nonzero limit, continue the first-order
radial modulation theory. If it vanishes, derive the distinguished radial
scale at which radial diffusion or higher radial derivatives enter before
building the finite-disk solver.

## Priority Gate R4: radial-diffusion promotion

Restore the cylindrical radial viscous and thermal diffusion terms in the
non-similar boundary-layer system and project their frozen-radius operators
onto the fold null pair. Compare

```text
chat_c_r, chat_c_rr, chat_c_r_geo
```

using the same wall-temperature transversality normalization. Separate the
scalar thermal Laplacian from the radial/azimuthal vector-Laplacian geometry
terms. Do not add axial-velocity diffusion to the continuity equation;
that requires a separate axial-momentum/pressure extension.

Decision: use the measured powers, not an assumed (c_{rr}=O(1)), to select
the radial scale. Then derive a variable-radius reduced amplitude equation
before committing to a global finite-disk solver.

## Priority Gate R5: two-dimensional thermo-radial outer

Use \(X=\epsilon R/\ell\) and \(Y=\epsilon z/\ell\) with the critical
fold--tail field scales. Derive the leading PDE directly from the
radially diffusive equations and require exact recovery of the
one-dimensional critical outer DAE when all \(X\)-derivatives vanish.
Separate the leading thermo-azimuthal system from the next meridional
momentum/pressure reconstruction.

Before solving the PDE, derive its wall/overlap, axis, radial-edge, and
far-field data. Independently select the amplitude, detuning, and forcing
scales; do not assume the fold quadratic term is leading merely because
the linear radial drift and diffusion balance.

## Priority Gate R6: full-equation and physical-realizability audit

Because \(R\sim z\sim\ell/\epsilon\), repeat the outer hierarchy from the
full axisymmetric continuity, three momentum equations, pressure and
energy equation. Include the first nonzero meridional momentum order and
its pressure-integrability condition even when it is smaller than the
leading centrifugal balance.

Express the physical window with
\(\Pi_R=\epsilon\sqrt{Re_R}\), the density contrast
\(|\beta\Delta T|\), confinement \(\Pi_H=\epsilon H_c/\ell\), and the
gravity-to-meridional ratio. Do not present the current
\(|\beta\Delta T|\simeq0.662\) critical state as quantitatively
Boussinesq-valid without a matched low-Mach variable-density calculation.

## Gate 3: BEK bifurcation map

Use pseudo-arclength continuation to compute self-consistent traditional and canonical base-flow branches as `Ro` and `Tw` vary. Construct the fold loci

```text
Tw_c^traditional(Ro),
Tw_c^canonical(Ro),
```

within the numerically converged and Boussinesq-relevant range. Search for fold-pair creation/annihilation and cusp points rather than assuming that the von Karman fold persists throughout the family.

Deliverables:

- converged branch and fold data in the `(Ro,Tw)` plane;
- fold/cusp detection and non-degeneracy checks;
- profiles and balance diagnostics on each relevant branch;
- explicit separation of self-consistent endpoint comparisons from frozen or term-toggle diagnostics.

Falsification check: if endpoint fold loci are indistinguishable throughout the valid range, revise the closure-induced-bifurcation hypothesis rather than extending to expensive stability sweeps.

## Gate 4: fold mechanism and adjoint structural sensitivity

At representative folds, compute right and left null vectors and predict the displacement caused by each density--acceleration contribution:

```text
dTw_c/deta_j = -<psi, partial_eta_j R>/<psi, partial_Tw R>.
```

Validate each derivative against small finite changes in `eta_j`. Use continuous `eta_j` paths only where needed to determine whether a fold shifts, exits the model-valid range or annihilates with another fold at a cusp.

Deliverables:

- term-by-term fold sensitivities versus `Ro`;
- finite-change linearity ranges and error tables;
- physical-space maps of the forcing terms and left/right null modes;
- an identified mechanism for any change in fold topology.

Stop condition: do not describe a single acceleration term as causal unless its adjoint prediction, finite-change response and endpoint comparison agree.

## Gate 5: matched low-Mach reference

For selected small-temperature-difference cases, compare the two Boussinesq endpoints with a self-consistent low-Mach compressible reference under matched geometry, rotation and property assumptions.

Deliverables:

- matched reference cases near and away from the predicted fold region;
- separation of density--acceleration closure, density-law, property and compressibility errors;
- evidence for or against describing the traditional fold as closure-induced.

Do not infer accuracy from the canonical derivation alone, and do not treat formal large-`Tw` continuation as quantitative Boussinesq validation.

## Gate 6: similarity-subspace temporal stability

Compute the temporal spectrum along representative branches and folds. Verify whether a single real eigenvalue crosses zero, whether stability exchanges at the folds, and how the recovery time changes with `Ro` and closure.

Deliverables:

- spectral-order convergence and constraint checks;
- leading eigenvalues and modes along the selected branches;
- a critical-slowing-down map;
- language explicitly restricted to axisymmetric similarity-subspace stability.

## Gate 7: local three-dimensional modal sensitivity

Only after the base-flow mechanism is established, add the same physically derived closure endpoints and diagnostic derivatives to the Blackburn/BEK perturbation operator. Implement the adjoint spatial eigenproblem and compare first-order eigenvalue shifts with direct endpoint or small-`eta_j` calculations for Type-I and Type-II modes.

Deliverables:

- direct/adjoint biorthogonality and original-pencil residual tests;
- term-resolved `delta alpha` predictions;
- selected neutral curves and critical-parameter validation across `Ro`;
- a two-mode or pseudospectral treatment near modal degeneracy if required;
- language restricted to local parallel-flow stability.

## Gate 8: observable error indicator

Test whether equation-level closure residuals combined with fold or modal amplification predict errors in observables such as `Tw_c`, `alpha`, `R_c` and `beta_c`.

Deliverables:

- a verified first-order observable-error estimate;
- a topology-warning criterion near folds, cusps or mode interaction;
- a map of where simplified closure results are insensitive or critically amplified;
- an honest statement of whether the indicator is empirical, asymptotic or a proved bound.

Do not present intermediate `eta_j` values as calibrated physical models. The practical comparison remains between physically derived endpoint closures.

## Gate 9: generality and manuscript scope

Apply the established mechanism to the rotor--stator system only after the BEK result is complete. Check whether confinement changes the fold mechanism without confusing similarity-subspace stability with full three-dimensional stability.

For a focused first manuscript, prioritise:

1. BEK derivation and endpoint models;
2. traditional/canonical fold and cusp map;
3. adjoint fold mechanism;
4. selected low-Mach validation.

The full Type-I/Type-II survey and rotor--stator generalisation may become separate studies if including them obscures the main closure-induced critical-amplification result.

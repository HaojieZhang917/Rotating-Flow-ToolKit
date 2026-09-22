# Corrected manual supplement calculation plan (2026-09-07)

## Purpose

Supply the minimum missing corrected-energy evidence needed to revise the learning manual after `LEARNING_MANUAL_AUDIT_2026-09-07.md`. The study tests finite-epsilon augmented-fold tangents across spatial degree. It does not compute a new two-parameter cusp curve or a three-dimensional stability problem.

## Model and parameters

- Acceleration-consistent steady BEK similarity model in `work/src/BEKConsistent.jl`.
- Energy equation `T'' + Pr*Ro*H*T' = 0`.
- `Pr=0.72`, `gamma=1`, rational map `(a,b,c)=(2,0.6,0.5)`.
- Degrees `N=200,240,280,320`.
- Physical fold coordinate `epsilon=-Hinf`; evaluate `epsilon=0.02,0.015,0.01`.
- Nonlinear tolerance `1e-9`; tangent is the implicit derivative of the complete augmented fold system with `dHinf/depsilon=-1`.

## Inputs and source dependencies

Saved corrected-energy fold profiles under `work/results/corrected_energy_epsilon_fold_chain/`. Each saved `epsilon=0.02` fold supplies the initial field and singular vector. The fold is corrected at each requested `epsilon`; no historical uncorrected-energy output enters the calculation.

## Outputs and claims

- Script: `work/scripts/corrected_manual_tangent_convergence_20260907.jl`.
- New output root: `work/results/corrected_manual_tangent_convergence_20260907/`.
- Per-degree tangent tables and scaled right-null-mode profiles.
- A cross-degree summary will quantify spatial spread and finite-epsilon distance of `dTw/dRo` from the analytic compatibility slope.

The results may support a finite-epsilon convergence window and component-amplitude scaling only where the cross-degree data agree. They will not by themselves establish the epsilon-to-zero absolute tangent coefficients or function-space profile convergence.

## Completed results

The calculation completed for every planned degree and epsilon. N=280 versus N=320 gives relative spreads of 0.154% at epsilon=0.02, 0.138% at 0.015, and 1.37% at 0.01 for both parameter tangent components. The first two points therefore have a useful finite-epsilon two-grid check; the smallest point does not support an absolute tangent limit.

At N=320, the direct/adjoint slopes are -1.0425152608, -1.0372555395 and -1.0320453221. Their relative distances from the endpoint compatibility slope -1.0217715606 are 2.030%, 1.515% and 1.005%, respectively. The direct and adjoint values agree within the same discrete system, while convergence to the endpoint slope remains finite-epsilon asymptotics.

The N=320 amplitude ratios `max|vH|/max|vG|` are 0.07183, 0.05381 and 0.03584; `max|vF|/max|vG|` are 0.001748, 0.001001 and 0.0004535. These support the component orders but do not establish profile convergence in a norm.

The complete evidence scope and conclusion were appended as L026 in `RESEARCH_PROGRESS_LOG.md` before use in the revised manual.

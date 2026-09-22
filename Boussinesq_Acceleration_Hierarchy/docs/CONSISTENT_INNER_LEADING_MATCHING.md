# Inner leading problem and its connection to the slow outer tail

Date: 2026-09-05.

## Scope and evidence dependencies

Derive the inner leading equations on eta=O(1) and their overlap conditions
with the L013 thermo-rotational outer system. No new BVP or continuation is
performed. Existing-data algebraic checks use only:

- `work/src/BEKConsistent.jl`: governing equations and wall conditions;
- `docs/CONSISTENT_OUTER_THERMAL_TAIL.md`: leading outer system;
- `work/results/consistent_fold_tail_scaling/polynomial_fits.csv`: the two
  quadratic limiting-coordinate estimates;
- `work/results/consistent_epsilon_fold_chain/convergence_N320/epsilon_fold_table.csv`
  and its saved profiles: fixed-eta trends along the existing five-point chain.

Diagnostics are implemented in Julia in
`work/scripts/check_consistent_inner_matching.jl`, with new output under
`work/results/consistent_inner_leading_matching/`. They test consistency of
the expansion, not existence or uniqueness of a matched solution. The prior
spectral-resolution scatter and extrapolation uncertainty remain applicable.

## Regular two-scale hypothesis

Let r=Ro, s=(2-r-r^2)/2, omega=s+r and chi(T)=1-gamma(T-1).
Assume r approaches a nonzero r_t and Tw approaches Tw_t with chi(Tw_t)>0.
On compact eta intervals use

    H=H0+o(1), F=F0+o(1), G=G0+o(1), T=T0+o(1).

In the outer layer use Y=epsilon eta and

    H=epsilon h(Y)+..., F=epsilon^2 f(Y)+...,
    G=g(Y)+..., T=Theta(Y)+... .

For the regular overlap considered here, the outer inlet limits h(0+),
g(0+), Theta(0+) and the needed derivatives are finite. The overlap means
epsilon->0 with 1<<eta<<1/epsilon, not eta->infinity at fixed epsilon.
The condition Hinf=-epsilon alone does not prove H=O(epsilon) everywhere:
the outer balance and regularity assumptions are essential.

## Inner equations

Primes in this section denote eta derivatives; all coefficients are at r_t.
The O(1) system is

    H0' + 2F0 = 0,
    F0'' + chi(T0){r_t(F0^2+H0 F0')-(s_t+r_t G0)^2/r_t}
          + omega_t^2/r_t = 0,
    G0'' + chi(T0){r_t(2F0 G0+H0 G0')+2s_t F0} = 0,
    T0'' + Pr Ro H0 T0' = 0.

Wall conditions are H0(0)=F0(0)=G0(0)=0 and T0(0)=Tw_t.
Regular overlap requires

    H0(infinity)=0, F0(infinity)=0,
    G0(infinity)=g_m, T0(infinity)=Theta_m,

with decaying inner gradients, where g_m=g(0+), Theta_m=Theta(0+).
The physical conditions G=1,T=1 belong to Y->infinity.

The zero-pumping inner condition also gives the necessary identity

    integral_0^infinity F0(eta) d eta = 0.

This identity permits both an identically zero radial field and a
sign-changing radial field; it does not prove uniqueness of the inner state.

## Temperature in the regular inner limit

Integrating the energy equation gives

    T0'(eta)=C0 exp(Pr integral_0^eta H0(x) dx).

For a regular localised inner velocity field with a finite limit of
integral_0^eta H0, bounded T0 and the overlap gradient T0'(infinity)=0
force C0=0. Thus

    T0(eta) = Tw_t,   Theta_m = Tw_t.

One way to obtain the same result from the finite-epsilon problem is

    T_eta(0)=(1-Tw)/I_epsilon,
    I_epsilon=integral_0^infinity exp(Pr integral_0^eta H_epsilon(x) dx) d eta.

If the regular outer layer has O(1) temperature amplitude and the integrated
inner axial flow is finite, its contribution to I_epsilon is O(1/epsilon).
Consequently T_eta(0)=O(epsilon), and T0 is constant on every fixed eta
interval. An O(1) wall-to-ambient temperature difference is then carried
predominantly by the long outer layer.

The integrability qualification matters. H0->0 by itself is insufficient:
H0~-a/eta may produce T0'~eta^(-Pr a), allowing a nonconstant bounded T0
when Pr a>1. Such a singular overlap is outside the regular h(0+) assumption
and must not be excluded by the far-field value of H alone.

## Reduced inner momentum BVP and overlap algebra

With chi_m=1-gamma(Tw_t-1), the inner velocity problem becomes

    H0'+2F0=0,
    F0''+chi_m{r_t(F0^2+H0 F0')-(s_t+r_t G0)^2/r_t}
           +omega_t^2/r_t=0,
    G0''+chi_m{r_t(2F0G0+H0G0')+2s_tF0}=0,

subject to the wall and zero-pumping overlap conditions above.
Taking the inner far-field radial balance yields

    chi_m(s_t+r_t g_m)^2=omega_t^2,
    g_m=[omega_t/sqrt(chi_m)-s_t]/r_t.

The root is selected continuously from the physical far field and presumes
chi remains positive. For gamma>0, Tw_t!=1 and omega_t!=0, g_m cannot equal 1.
Theta_m and g_m are therefore not two independent matching constants:
Theta_m=Tw_t and g_m follows from radial balance, under the regular hypothesis.

The corresponding leading outer inlet data are Theta(0+)=Tw_t,
g(0+)=g_m, with the outer physical limits Theta(infinity)=1,g(infinity)=1,
h(infinity)=-1. Matching H0(infinity)=0 does not determine h(0+):
that coefficient belongs to the first nonzero axial-flow order. In particular,
the wall condition H(0)=0 cannot simply be imposed as h(0+)=0 without analysing
the inner correction.

## Possible degeneration of the O(1) inner core

If chi_m s_t^2=omega_t^2 (s_t!=0), the inner system admits the exact state

    F0=H0=G0=0, T0=Tw_t, g_m=0,
    Tw_t=1+[1-(omega_t/s_t)^2]/gamma.

For gamma=1 this is Tw_t=2-(omega_t/s_t)^2. It is a compatibility curve,
not an equation selecting r_t or a proof that a fold occurs there.
G0=0 denotes co-rotation with the disk in the adopted variables; it does
not mean that the dimensional absolute angular velocity vanishes.

If this state is selected, the leading inner momentum field is trivial and
the O(1) azimuthal change from 0 to 1 is carried by the outer layer.
Nonzero inner corrections can still be required for derivative and axial-flow
matching. Existing-coordinate and fixed-eta checks below assess this possibility;
they cannot prove that the full branch selects it or that no other inner
solutions exist.

## Existing-data check and evidence scope

The Julia diagnostic gives:

| limiting fit | inferred g_m | Tw_t minus zero-core compatibility temperature |
|---|---:|---:|
| five-point quadratic | -2.9791e-5 | 8.7059e-6 |
| four-point quadratic | -3.3629e-6 | 9.8275e-7 |

These defects are smaller than the change of fitted intercepts between the
two windows. The small negative g_m is not evidence for a resolved reversed
inner plateau; it is consistent with g_m=0 within the extrapolation uncertainty.

At fixed eta=1 the saved N=320 fields give:

| epsilon | H | F | G | Tw-T |
|---:|---:|---:|---:|---:|
| 0.100 | -1.0167e-3 | 1.4839e-3 | 0.074693 | 0.024791 |
| 0.050 | -2.8529e-4 | 3.8910e-4 | 0.038445 | 0.011939 |
| 0.020 | -4.7809e-5 | 6.3330e-5 | 0.015615 | 0.0046729 |
| 0.010 | -1.2104e-5 | 1.5898e-5 | 0.0078452 | 0.0023198 |
| 0.005 | -3.0820e-6 | 3.9309e-6 | 0.0039318 | 0.0011558 |

The fixed-eta trend supports vanishing leading inner velocities and a constant
inner temperature. It does not identify an overlap plateau at finite epsilon
or replace a uniform convergence proof. It suggests G and the inner temperature
drop are first-order corrections, while H and F at this station are smaller;
no new correction exponent is established here.

Raw diagnostics: `work/results/consistent_inner_leading_matching/limit_compatibility.csv`
and `fixed_eta_profiles.csv`. Source interpolations are piecewise linear and
all conclusions inherit the previous numerical-resolution limits. See log L014.

## What the leading matching establishes

Under the stated regularity hypothesis:

    (H0,F0,G0,T0)(infinity)=(0,0,g_m,Tw_t),
    (h,f,g,Theta)(infinity)=(-1,0,1,1),
    g(0+)=g_m, Theta(0+)=Tw_t.

It leaves h(0+), flux amplitudes and higher-order inner corrections to the
next solvability calculation. The slow manifold is an algebraic momentum
constraint here; no normal-hyperbolicity or dynamical-attraction theorem
has been proved by calling it a manifold.

The appropriate leading picture can therefore contain a trivial O(1) inner
core and a nontrivial O(1/epsilon) thermo-rotational tail. The overlap is an
asymptotic common range, not a third independently specified physical layer.

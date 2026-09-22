# Radial transverse response of heated-BEK similarity solutions

## Purpose and terminology

The immediate theoretical question is whether the heated-BEK similarity
invariant subspace is radially attracting, repelling, or singularly
sensitive to non-similar forcing. This is a steady spatial-response problem,
not temporal stability and not full three-dimensional stability.

Use \(\xi=\log x\), so \(x\partial_x=\partial_\xi\), and perturb a
similarity solution \(Q_0(\eta)=(H_0,F_0,G_0,T_0)\) by

\[
Q(\xi,\eta)=Q_0(\eta)+\delta e^{\lambda\xi}
(h,f,g,\theta)(\eta)+O(\delta^2).
\]

For increasing radius, \(\Re\lambda<0\) denotes outward decay and
\(\Re\lambda>0\) outward growth. This is a spatial exponent, not a temporal
growth rate.

## Continuous transverse pencil

Define

\[
q_0=s+RoG_0,\quad \chi_0=1-\gamma(T_0-1),
\]

\[
A_{r0}=Ro(F_0^2+H_0F_0')-{q_0^2\over Ro},
\quad
A_{\theta0}=Ro(2F_0G_0+H_0G_0')+CoF_0.
\]

Linearization of the non-similar PDE gives

\[
\boxed{
\begin{aligned}
h'+(2+\lambda)f&=0,\\
f''+\chi_0\left[
Ro\{(2+\lambda)F_0f+H_0f'+F_0'h\}-2q_0g
\right]-\gamma A_{r0}\theta&=0,\\
g''+\chi_0\left[
Ro\{2G_0f+(2+\lambda)F_0g+H_0g'+G_0'h\}
+Co f\right]-\gamma A_{\theta0}\theta&=0,\\
\theta''+Pr\,Ro\left(H_0\theta'+T_0'h+\lambda F_0\theta\right)&=0.
\end{aligned}}
\]

For uniform wall temperature, use

\[
h(0)=f(0)=g(0)=\theta(0)=0,
\qquad
f(\infty)=g(\infty)=\theta(\infty)=0.
\]

No independent condition is imposed on \(h(\infty)\); it is the perturbation
of the emergent axial flux. A fixed-flux apparatus closure would define a
different spatial spectrum.

In operator form,

\[
\boxed{\left(J_{\rm sim}+\lambda M_r\right)\widehat Q=0.}
\]

In the field order \((h,f,g,\theta)\), the only nonzero interior blocks of
\(M_r\) are

\[
(M_r)_{H,F}=I,\qquad
(M_r)_{F,F}=\chi_0RoF_0,
\]

\[
(M_r)_{G,G}=\chi_0RoF_0,\qquad
(M_r)_{T,T}=Pr\,RoF_0.
\]

All physical boundary rows of \(M_r\) vanish.

## Connection to the similarity fold

At \(\lambda=0\), the pencil becomes

\[
J_{\rm sim}\widehat Q=0.
\]

Thus a simple fixed-\(T_w\) similarity fold supplies a transverse spatial
mode at \(\lambda=0\). Away from a fold, an exponent approaching zero
signals increasingly long radial response. Relevant observables are:

1. finite generalized eigenvalues of \((-J_{\rm sim},M_r)\);
2. the exponent closest to zero and its right/left mode;
3. the resolvent norm for imposed radial wall heating;
4. their behavior approaching a similarity fold or thermal tail.

The identity does not prove transverse stability of the full PDE. It proves
that similarity folds are zero-spatial-exponent singularities of this
parabolized radial-response model.

## Fold-centered radial amplitude equation

The local geometry can be exposed without first computing a complete
finite-disk field. Let \(U_f\) be a simple similarity fold and let

\[
J_fv=0,\qquad J_f^\dagger w=0.
\]

Write a slowly modulated non-similar state as

\[
U(\xi,\eta)=U_f(\eta)+A(\xi)v(\eta)+\cdots,
\]

and allow

\[
\Delta T_w=T_w-T_{w,f},\qquad
\Delta Ro=Ro-Ro_f.
\]

Expanding the steady radial PDE through the first fold-normal-form order
and projecting with \(w\) gives

\[
\boxed{
c_r A_\xi+c_T\Delta T_w+c_R\Delta Ro+c_2A^2
=\mathcal F_{\rm ns}(\xi)+\cdots,}
\]

where

\[
\begin{aligned}
c_r&=\langle w,M_rv\rangle,\\
c_T&=\langle w,R_{T_w}\rangle,\\
c_R&=\langle w,R_{Ro}\rangle,\\
c_2&={1\over2}\langle w,R_{UU}[v,v]\rangle .
\end{aligned}
\]

The forcing \(\mathcal F_{\rm ns}\) represents a specified radial wall
temperature, edge, or upstream perturbation. The same continuous pairing
and boundary contributions must be used in all four coefficients.

For \(A_\xi=0\), this reduces to the ordinary similarity fold normal form.
Thus \(c_2\) measures the local curvature of the similarity branch, while
\(c_r\) converts that fold geometry into radial response. On a nearby
equilibrium branch \(A=A_*\),

\[
\lambda_{\rm slow}=-{2c_2A_*\over c_r}+\cdots.
\]

This is the small spatial exponent expected from the full pencil. At the
fold itself,

\[
A_\xi=-{c_2\over c_r}A^2+\cdots,
\]

so the response is algebraic in \(\xi=\log x\), rather than exponentially
attracting or repelling.

Three outcomes have distinct meanings:

1. \(c_r\neq0\): the parabolized PDE has a first-order radial fold
   unfolding, and the sign of \(2c_2A_*/c_r\) determines the local radial
   attraction direction;
2. \(c_r=0\) but a higher radial coefficient is nonzero: the first-order
   transverse problem is degenerate and radial diffusion or a higher-order
   modulation equation is leading;
3. \(c_r\) tends to zero in the fold--tail limit: radial non-similarity and
   tail delocalization interact in a new distinguished limit.

This amplitude equation is the precise mathematical version of asking
whether the similarity solution lies in a locally attracting or repelling
radial channel. Linear response determines the local channel; it does not
exclude disconnected non-similar steady branches.

## Axis and edge interpretation

Analytic axisymmetric perturbations have even radial expansions, so regular
single-axis contributions have \(\lambda=2n\), \(n=0,1,\ldots\). General
complex \(\lambda\) modes are naturally Mellin modes on an annulus away from
the axis. Edge responses require their superposition and eventually radial
diffusion, which is absent from this reduced pencil.

This theory therefore measures whether the similarity interior is weakly
or strongly conditioned to long-wave radial forcing. It does not replace
the global finite-disk calculation.

## Required numerical checks

Before interpreting any \(\lambda\):

1. compare \(J_{\rm sim}+\lambda M_r\) with a centered directional
   derivative of the nonlinear non-similar residual;
2. separate finite generalized eigenvalues from infinite eigenvalues caused
   by singular \(M_r\) and boundary rows;
3. vary the wall-normal grid and rational map independently;
4. check direct and adjoint pencil residuals;
5. distinguish fixed-\(T_w\), forced-\(T_w(x)\), and fixed-flux closures;
6. never call \(\Re\lambda\) a temporal growth rate.

The first study will use \(Ro=-0.50,-0.4703674,-0.40\), but each state must
also have an explicitly documented \(T_w\) or \(\epsilon\). Rossby number
alone does not identify a nonlinear branch member.

## Pencil verification result

The shared Julia assembly is in
'work/src/BEKRadialResponse.jl'. An independent nonlinear residual check in
'work/scripts/verify_radial_response_pencil_20260909.jl' perturbs a
similarity base by a radial Mellin mode and compares a centered directional
difference with \((J_{\rm sim}+\lambda M_r)\widehat Q\).

Using degree 60 isothermal consistent bases at the three provisional
Rossby numbers and

\[
\lambda=-1.3,\ -0.2,\ 0,\ 0.4,\ 2,
\]

the relative infinity-norm discrepancies range from
\(2.23\times10^{-12}\) to \(6.55\times10^{-11}\). This passes the
operator-assembly Gate, including the \(\lambda=0\) similarity-Jacobian
limit. These are derivative-verification data, not transverse-spectrum
results; the isothermal bases were chosen only because they are regular
regression states.

The recorded output is
'work/results/nonsimilar_heated_bek_20260909/radial_pencil_directional_check.csv'.

## Fold-chain radial-transversality Gate

The fold-chain calculation must not compare the raw value obtained after
imposing $w^Tv=1$. On the saved collocation folds, $w^Tv$ becomes small
and changes sign at some resolutions, so that convention can create a pole
in both $c_r$ and $c_2$. A normalization-invariant diagnostic is

\[
 \widehat c_r={w^TM_rv\over w^TR_{T_w}},\qquad
 \widehat c_2={\tfrac12w^TR_{UU}[v,v]\over w^TR_{T_w}}.
\]

The right mode is fixed by $\max_\eta|v_G|=1$, with a fixed sign. Scaling
$w$ cancels from both ratios, and the normal form becomes

\[
 \widehat c_r A_\xi+\Delta T_w+\widehat c_2A^2=\cdots.
\]

The corrected-energy folds give:

| ϵ | ĉ_r, N=240 | ĉ_r, N=280 | ĉ_r, N=320 |
|---:|---:|---:|---:|
| 0.10 | -0.0244691 | -0.0244664 | -0.0244672 |
| 0.05 | -0.0117985 | -0.0117992 | -0.0117996 |
| 0.02 | -0.00459642 | -0.00459887 | -0.00459642 |
| 0.01 | -0.00227710 | -0.00227738 | -0.00227773 |
| 0.005 | -0.00114324 | -0.00113150 | -0.00113160 |

Here ĉ denotes the normalization-invariant hatted coefficient. A
log--log fit over ϵ=0.10,0.05,0.02,0.01 gives, at both $N=280$ and
$320$,

\[
 \boxed{|\widehat c_r|\simeq0.260\,\epsilon^{1.030}}.
\]

Thus the present evidence answers the radial-transversality Gate in the
affirmative: relative to the ordinary fold unfolding, the first-order
radial unfolding vanishes linearly as the fold approaches the tail limit.

For the same four resolved points, the $N=320$ provisional fits are

\[
 |\widehat c_2|\simeq3.96\,\epsilon^{2.145},\qquad
 |c_2/c_r|\simeq15.2\,\epsilon^{1.114}.
\]

The $c_2$ value at ϵ=0.005 is not converged and is excluded, so the
robust conclusion concerns $\widehat c_r=O(\epsilon)$, not yet the precise
power of $c_2$.

### Direct pencil check away from the fold

For $N=200$, two nearby fixed-$T_w$ branches were computed at each of
ϵ=0.10,0.05,0.02. A bordered eigenpair corrector solved
$(J+\lambda M_r)q=0$ from the fold mode. The direct slow spatial
eigenvalue differs from

\[
 \lambda_{\rm NF}=-{2c_2A_*\over c_r}
\]

by 1.2--3.0 percent across all six tests and changes sign between the two
branches. This independently validates the local coefficient combination.

The calculation is in
`work/scripts/analyze_fold_radial_transversality_20260909.jl`,
`work/scripts/validate_fold_radial_normal_form_20260909.jl`, and
`work/results/fold_radial_transversality_20260909/`.

The current comparison varies $N$ at the established rational map
$(a,b,c)=(2,0.6,0.5)$. Mapping sensitivity remains outstanding, so the
precise prefactor is not yet a final continuum value.

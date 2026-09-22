# Non-similar heated BEK boundary-layer model (2026-09-09)

## Research question and evidence status

The purpose of this model is to test

\[
\textit{Does the fold--tail topology of the heated BEK similarity equations
survive when radial thermal development is restored?}
\]

This document defines and checks the model before any numerical result is
used. No two-dimensional solution, fold, or topology claim is made here.
The previous matched-adjoint calculation is paused until this physical
validation gate is resolved.

The first benchmark set is provisionally

\[
Ro=-0.50,\qquad -0.4703674,\qquad -0.40,
\quad Pr=0.72,\quad \gamma=1,
\]

but these numbers are not sufficient to define a non-similar problem:
radial wall and inflow/edge data must also be specified. New numerical
output will be written under
'work/results/nonsimilar_heated_bek_20260909/'; nothing in 'baselines/'
will be modified.

## Dimensional parent boundary-layer equations

Use cylindrical inertial coordinates \((R,\phi,z)\), assume steady
axisymmetry, and write

\[
u_R=R\,A(R,z),\qquad u_\phi=R\,B(R,z),\qquad u_z=W(R,z),
\qquad T=T(R,z).
\]

The temperature is a scalar and is not multiplied by \(R\). Let
\(\chi(T)=\rho/\rho_\infty=1-\gamma(T-1)\). In the acceleration-consistent
Boussinesq closure, \(\chi\) multiplies the complete material acceleration,
while the imposed outer pressure gradient remains the reference-density
gradient. With the boundary-layer approximation retaining only
wall-normal diffusion,

\[
\begin{aligned}
2A+R A_R+W_z&=0,\\
\chi\left(A^2+RAA_R+WA_z-B^2\right)
&=-\Omega_F^2+\nu A_{zz},\\
\chi\left(2AB+RA B_R+WB_z\right)&=\nu B_{zz},\\
RA T_R+WT_z&=\alpha T_{zz}.
\end{aligned}
\]

The radial pressure relation used here is
\(p_R/\rho_\infty=\Omega_F^2R\), inherited from the rigidly rotating
far field. Radial viscous and thermal diffusion, streamwise pressure
adjustment, and axial pressure variation are omitted in this first
parabolized model. These are model assumptions, not consequences of
similarity.

## BEK nondimensional form

Let

\[
\ell=\sqrt{\nu/\Omega},\qquad x={R\over L},\qquad \eta={z\over\ell},
\]

\[
s={\Omega_D\over\Omega}={Co\over2},\qquad
Ro={\Omega_F-\Omega_D\over\Omega},\qquad
\omega=s+Ro={\Omega_F\over\Omega},
\]

and introduce

\[
u_R=-\Omega Ro\,R\,F(x,\eta),\qquad
u_\phi=\Omega R\{s+Ro\,G(x,\eta)\},
\qquad
u_z=-\ell\Omega Ro\,H(x,\eta).
\]

Here \(L\) only defines the radial coordinate; because radial diffusion is
discarded, it disappears after \(R\partial_R=x\partial_x\). With
\(q=s+RoG\) and \(\chi=1-\gamma(T-1)\), the equations are

\[
\boxed{
\begin{aligned}
H_\eta+2F+xF_x&=0,\\
F_{\eta\eta}
+\chi\left[
Ro\{F^2+xFF_x+HF_\eta\}-{q^2\over Ro}
\right]+{\omega^2\over Ro}&=0,\\
G_{\eta\eta}
+\chi\left[
Ro\{2FG+xFG_x+HG_\eta\}+CoF
\right]&=0,\\
T_{\eta\eta}
+Pr\,Ro\{xFT_x+HT_\eta\}&=0.
\end{aligned}}
\]

Thus radial heat advection is restored as \(Pr\,Ro\,xFT_x\). Its
dimensional magnitude relative to axial heat advection is

\[
\mathcal R_T(x,\eta)={|xF T_x|\over |H T_\eta|},
\]

where a signed budget or regularized norm must replace the ratio near zeros
of the denominator.

Candidate boundary data are

\[
F=0,\quad G=0,\quad H=0,\quad T=T_w(x)\qquad(\eta=0),
\]

\[
F\to0,\quad G\to1,\quad T\to1\qquad(\eta\to\infty).
\]

A radial inflow profile is additionally required wherever the system is
marched into the domain.

## Gate 1: exact recovery of the similarity equations

Set \(F_x=G_x=H_x=T_x=0\). The PDE reduces exactly to

\[
\begin{aligned}
H'+2F&=0,\\
F''+\chi\left[Ro(F^2+HF')-{(s+RoG)^2\over Ro}\right]
+{\omega^2\over Ro}&=0,\\
G''+\chi\left[Ro(2FG+HG')+CoF\right]&=0,\\
T''+Pr\,RoHT'&=0,
\end{aligned}
\]

which is the residual implemented in 'work/src/BEKConsistent.jl'.
Consequently,

\[
\boxed{\text{heated BEK similarity solutions form an exact invariant
subspace of the proposed PDE.}}
\]

This does not establish whether that subspace attracts, repels, or
approximates solutions with non-similar radial data.

## Gate 2: radial marching direction and degeneracy

The dimensional radial transport velocity is

\[
u_R=-\Omega Ro\,R F.
\]

Therefore marching direction is controlled by \(-RoF\), not by \(F\)
alone. For \(Ro<0\), \(u_R\) has the same sign as \(F\). If this sign is
uniform through the boundary layer, information propagates uniformly
outward or inward. If it changes with \(\eta\), a one-way radial march
does not supply all characteristics and a global two-dimensional
boundary-value formulation may be required.

The logarithmic coordinate \(\xi=\log x\) gives
\(x\partial_x=\partial_\xi\). The coefficients of
\(F_\xi,G_\xi,T_\xi\) contain \(F\), so \(F=0\) is a characteristic
degeneracy. It occurs at the wall and far field even if \(F\) has one sign
in the interior. A method based on pointwise division by \(F\) is therefore
not admissible.

Axis regularity gives

\[
\begin{aligned}
F&=F_0(\eta)+x^2F_2(\eta)+\cdots,&
G&=G_0(\eta)+x^2G_2(\eta)+\cdots,\\
H&=H_0(\eta)+x^2H_2(\eta)+\cdots,&
T&=T_0(\eta)+x^2T_2(\eta)+\cdots.
\end{aligned}
\]

At \(x=0\), all \(x\partial_x\) terms vanish and the leading profiles
satisfy the similarity ODE. The \(O(x^2)\) system supplies an axis launch
only after radial wall or outer data determine its forcing.

## A necessary non-similarity condition

For an infinite disk with constant \(T_w\), rigid wall and far-field
rotation, and radially uniform outer conditions, all boundary data admit
\(F_x=G_x=H_x=T_x=0\). Starting from an exact similarity profile therefore
reproduces that solution and gives \(\mathcal R_T=0\) identically. Such a
calculation is only a code regression test; it cannot answer whether radial
development destroys the fold or tail.

A discriminating benchmark must specify at least one source of radial
development:

1. a finite-radius disk with an edge and radial inflow data;
2. a prescribed nonuniform wall temperature \(T_w(x)\);
3. non-similar upstream profiles corresponding to a stated apparatus.

These are different physical problems and cannot be inferred from \(Ro\)
and one scalar \(T_w\). Until one is fixed, the proposed three-\(Ro\)
sweep is under-specified.

## Minimal numerical gates

The solver will be written in Julia under 'work/src/' and its driver under
'work/scripts/'. Before a topology sweep it must pass:

1. constant radial data reproduce the existing ODE profiles and all radial
   advection budgets vanish;
2. a manufactured solution verifies the signs of \(xF_x,xFG_x,xFT_x\);
3. the sign of \(-RoF(x,\eta)\) and its zero contours are recorded;
4. the \(\eta\) grid, radial step, start radius, and outer truncation are
   varied independently;
5. signed and absolute energy contributions from \(xFT_x\), \(HT_\eta\),
   and \(T_{\eta\eta}/(PrRo)\) are reported.

Only after these gates pass should the three provisional \(Ro\) values be
heated and compared with the similarity fold/tail topology.

## Completed Gate results

The Julia similarity check evaluated the PDE and ODE residuals at 1000
random local states for each provisional Rossby number, with all radial
derivatives set to zero. The largest difference was
\(4.441\times10^{-16}\), which is floating-point roundoff. Gate 1 passes.

The direction audit used the existing no-fold first-order generic-tail
inner profiles. On \(0<\eta\leq8\), the coefficient
\(u_R/(\epsilon\Omega R)=-RoF_1+O(\epsilon)\) behaves as follows:

| \(Ro\) | minimum | maximum | interpretation |
|---:|---:|---:|---|
| \(-0.50\) | \(-2.4022\times10^{-2}\) | \(1.0381\times10^{-3}\) | sign change near \(\eta=5.3315\) |
| \(-0.470367\) | \(-5.8\times10^{-12}\) | \(1.23\times10^{-10}\) | first-order field is degenerate |
| \(-0.40\) | \(-1.6710\times10^{-3}\) | \(3.8669\times10^{-2}\) | sign change near \(\eta=4.6177\) |

The near-\(Ro_t\) signs are numerical noise around a vanishing first-order
field. On both generic sides, different parts of the same radial station
receive information from opposite radial directions. A scalar one-way
march cannot therefore be the primary three-case validation solver.

## Consequence for the finite-disk model

A finite disk remains the preferred physical symmetry breaking, but its
edge cannot generally be represented by an imposed downstream condition in
the parabolized equations. Near the edge, radial gradients have the
boundary-layer scale and the discarded radial diffusion terms return. For
the angular-rate fields and temperature, they are

\[
\delta_R^2\left(F_{xx}+{3\over x}F_x\right),\qquad
\delta_R^2\left(G_{xx}+{3\over x}G_x\right),\qquad
\delta_R^2\left(T_{xx}+{1\over x}T_x\right),
\quad \delta_R={\ell\over L}.
\]

Thus an edge-invasion calculation requires a global axisymmetric model
retaining these terms, or the existing full compressible solver with
matched geometry and boundary data. What lies outside the disk must also be
represented: the switch from a rotating heated no-slip surface to the
exterior is an edge-region problem, not one extra marching condition.

The revised hierarchy is:

1. retain the parabolized PDE as an algebraic regression and controlled
   weak-\(T_w(x)\) experiment;
2. use a global finite-disk model for the primary physical test;
3. vary \(R_d/\ell\) and \(\epsilon=-H_\infty\) independently to test
   whether the large-radius and tail limits commute;
4. resume the matched-adjoint calculation only if the global solution has
   an interior region approaching the similarity fold--tail structure.

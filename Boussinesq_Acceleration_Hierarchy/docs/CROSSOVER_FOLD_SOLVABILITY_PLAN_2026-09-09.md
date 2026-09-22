# Crossover fold-solvability plan (2026-09-09)

## Objective

The basic-flow calculation determines the matching-degeneracy location
$Ro_t$. The remaining problem is to determine the local fold direction

\[
\Lambda_f=\lim_{\epsilon\to0}{Ro_f(\epsilon)-Ro_t\over\epsilon}
\]

from the fold equations rather than from a coordinate fit.

## Actual augmented system

The code solves, at prescribed $H_\infty=-\epsilon$,

\[
R(U,T_w,Ro)=0,\qquad
J(U,T_w,Ro)v=0,\qquad
\langle v,v\rangle=1,
\]

where $J=R_U$. The equation fixing $H_\infty$ releases $Ro$. Thus the
fold condition is not an extra condition on the leading basic-flow outer
problem; it enters through the null equation and its next-order
compatibility.

In the crossover layer write

\[
Ro=Ro_t+\epsilon\Lambda+\cdots,
\qquad
U=U_0+\epsilon U_1(\Lambda)+\epsilon^2U_2(\Lambda)+\cdots,
\]

and expand the nullvector and Jacobian as

\[
v=v_0+\epsilon v_1+\cdots,
\qquad
J=\mathcal L_0+\epsilon
\left(\mathcal L_1^{(0)}+\Lambda\mathcal L_1^{(\mu)}\right)+\cdots.
\]

The leading and first-correction null equations are

\[
\mathcal L_0v_0=0,
\]

and

\[
\mathcal L_0v_1+
\left(\mathcal L_1^{(0)}+\Lambda\mathcal L_1^{(\mu)}\right)v_0=0.
\]

If $w_0$ spans the adjoint nullspace, the Fredholm condition is

\[
\boxed{
\mathcal S(\Lambda)=
\langle w_0,\mathcal L_1^{(0)}v_0\rangle+
\Lambda\langle w_0,\mathcal L_1^{(\mu)}v_0\rangle=0,
}
\]

so that, for a simple transverse root,

\[
\boxed{
\Lambda_f=-
{\langle w_0,\mathcal L_1^{(0)}v_0\rangle
 \over
 \langle w_0,\mathcal L_1^{(\mu)}v_0\rangle}.}
\]

Because the limiting operator is an inner--outer matched operator, these
inner products cannot be replaced without proof by an unweighted Euclidean
dot product of finite collocation vectors.

## Required preliminary null-mode scaling check

The crossover base state has $F,H=O(\epsilon^2)$ and $G,T-T_t=O(\epsilon)$.
The existing fold modes suggest the distinct null scaling

\[
v_G,v_T=O(1),\qquad v_H=O(\epsilon),\qquad v_F=O(\epsilon^2).
\]

Before deriving the limiting adjoint problem, this will be checked using
the corrected full fold modes at $N=200,240,280,320$ and
$\epsilon=0.02,0.015,0.01$. The source data are
`work/results/corrected_manual_tangent_convergence_20260907/` and
`work/results/presentation_evidence_audit_20260908/tangent/`.
Post-processed output will be written to
`work/results/crossover_fold_solvability_20260909/`.

This diagnostic can validate the form of the expansion, but it cannot by
itself evaluate the matched Fredholm quotient above.

## Preliminary null-mode result

The corrected fold modes at $N=240,280,320$ collapse under the proposed
componentwise scaling. For $N=320$, linear extrapolation using
$\epsilon=0.02,0.015,0.01$ gives

\[
{\max|v_H|\over\epsilon\max|v_G|}\longrightarrow3.5766,
\qquad
{\max|v_F|\over\epsilon^2\max|v_G|}\longrightarrow4.6993,
\]

and

\[
{\max|v_T|\over\max|v_G|}\longrightarrow0.7963.
\]

The three grids agree closely for the scaled amplitudes; the $N=240$
$dRo/d\epsilon$ value at $\epsilon=0.01$ is an outlier, while $N=280$ and
$N=320$ remain close. A direct linear fit of the three N320 tangents gives
$dRo/d\epsilon\to-0.3070$ with RMS 0.0021. This is consistent with a finite
fold direction near $-0.30$, but it is not yet a reliable determination of
$\Lambda_f$: coordinate-quotient and derivative extrapolations retain
finite-$\epsilon$ and grid sensitivity.

The key result here is therefore the scaling, not the extrapolated decimal:

\[
\boxed{v_G,v_T=O(1),\qquad v_H=O(\epsilon),\qquad v_F=O(\epsilon^2).}
\]

This supplies a nontrivial limiting $v_0$ for the next-order Fredholm
problem. The next calculation must construct the consistently scaled left
null mode and determine how its mass separates between the inner and outer
regions before evaluating $\mathcal S_0$ and $\mathcal S_1$.

## Continuous outer limiting null operator

Let $(h,f,g,\theta)$ denote the critical leading outer solution, with

\[
q=s_t+Ro_tg,\qquad
\chi=1-\gamma(\theta-1),\qquad
B=Ro_t(2fg+hg')+2s_tf.
\]

Use outer null variables $(a,b,c,d)$ defined by

\[
v_H=\epsilon a(Y),\quad
v_F=\epsilon^2b(Y),\quad
v_G=c(Y),\quad
v_T=d(Y),\qquad Y=\epsilon\eta.
\]

Linearizing the continuous outer DAE at fixed $Ro_t$ gives

\[
\begin{aligned}
\mathscr C&=a'+2b,\\
\mathscr A&=2\chi qRo_t c-\gamma q^2d,\\
\mathscr G&=c''+\chi Ro_thc'+2\chi Ro_tf c
 +\chi Ro_tg'a+2\chi q b-\gamma B d,\\
\mathscr T&=d''+\mathrm{Pr}\,Ro_thd'+\mathrm{Pr}\,Ro_t\theta'a.
\end{aligned}
\]

Thus the outer limiting null problem is a differential--algebraic problem,
not four uncoupled differential equations. In particular, the radial slow-
manifold constraint fixes the leading relation between $c$ and $d$.

## Outer Green identity and formal adjoint

Introduce multipliers $(\ell,\rho,\sigma,\tau)$ for
$(\mathscr C,\mathscr A,\mathscr G,\mathscr T)$. Direct integration by parts
gives

\[
\int\left(\ell\mathscr C+\rho\mathscr A+\sigma\mathscr G
+\tau\mathscr T\right)dY
=\int\left(a\mathscr C^\dagger+b\mathscr A^\dagger
+c\mathscr G^\dagger+d\mathscr T^\dagger\right)dY
+[\mathscr B_o],
\]

where

\[
\begin{aligned}
\mathscr C^\dagger
 &=-\ell'+\chi Ro_tg'\sigma+\mathrm{Pr}\,Ro_t\theta'\tau,\\
\mathscr A^\dagger
 &=2\ell+2\chi q\sigma,\\
\mathscr G^\dagger
 &=2\chi qRo_t\rho+\sigma''-(\chi Ro_th\sigma)'
   +2\chi Ro_tf\sigma,\\
\mathscr T^\dagger
 &=-\gamma q^2\rho-\gamma B\sigma+\tau''
   -(\mathrm{Pr}\,Ro_th\tau)'.
\end{aligned}
\]

The outer boundary concomitant is

\[
\boxed{
\mathscr B_o=
\ell a+\sigma c'-\sigma'c+\chi Ro_th\sigma c
+\tau d'-\tau'd+\mathrm{Pr}\,Ro_th\tau d.}
\]

At the far field, the primal conditions make $a,c',d'$ free while
$c=d=0$; hence a separated outer adjoint would require
$\ell=\sigma=\tau=0$ there. At $Y\downarrow0$, however, $c'$ and $d'$ are
matching data rather than arbitrary boundary data. The entrance value of
$\mathscr B_o$ must therefore be cancelled by the inner overlap
concomitant, not set to zero prematurely.

## Why one global variable rescaling is insufficient

The max-norm scaling above is a global statement dominated by the outer
tail. Since the outer null fields satisfy zero entrance values, their inner
limits have an additional power of $\epsilon$. At fixed $\eta$ the
consistent inner variables are instead

\[
v_G=\epsilon\bar c(\eta),\qquad
v_T=\epsilon\bar d(\eta),\qquad
v_H=\epsilon^2\bar a(\eta),\qquad
v_F=\epsilon^2\bar b(\eta).
\]

Their lowest equations are

\[
\bar a'+2\bar b=0,\qquad
-2\chi_ts_tRo_t\bar c+\gamma s_t^2\bar d=0,\qquad
\bar c''=0,\qquad
\bar d''=0.
\]

If these equations alone are adjointed, the multiplier of the continuity
equation is forced to zero and the radial viscous term for $\bar b$ is
absent. That is not the full matched Fredholm operator: the missing radial
equation appears one order later and is precisely of the order needed by
the crossover correction. Therefore both variable columns and residual
rows must be asymptotically weighted. The correct object is a graded
inner--outer operator, schematically

\[
\widetilde{\mathcal L}_\epsilon
=D_R(\epsilon)^{-1}\mathcal L_\epsilon D_v(\epsilon),
\]

with different inner and outer versions joined by matching. This prevents
the leading inner limit from discarding the equation that carries the
solvability information.

The next derivation must retain that next radial inner row, derive its
concomitant, and show explicitly that its overlap contribution cancels the
$Y\downarrow0$ value of $\mathscr B_o$.

The regional scaling is also directly visible in the corrected fold-mode
profiles. For $N=320$ at fixed $\eta=1$, as
$\epsilon:0.02\to0.01$,

\[
{v_G\over\epsilon\max|v_G|}:2.0267\to2.0403,\qquad
{v_T\over\epsilon\max|v_G|}:-0.5883\to-0.5848,
\]

and

\[
{v_H\over\epsilon^2\max|v_G|}:-2.5556\to-2.5829,\qquad
{v_F\over\epsilon^2\max|v_G|}:2.4329\to2.4561.
\]

At fixed $Y=0.5$, the outer-scaled quantities instead converge directly;
for example $v_H/(\epsilon\max|v_G|)$ changes only from $-2.6869$ to
$-2.6891$. This confirms from the same global eigenfunction that the two
regions require different component scalings.

## Graded inner null hierarchy through the radial-viscous order

Write the crossover base state in fixed inner coordinates as

\[
F=\epsilon^2f_2+\cdots,\quad H=\epsilon^2h_2+\cdots,\quad
G=\epsilon g_1+\cdots,\quad
T=T_t+\epsilon\theta_1+\cdots,
\]

where

\[
g_1=p_t\eta,\qquad
\theta_1=\mathcal T'(Ro_t)\Lambda+B_T\eta.
\]

For the null mode use

\[
v_G=\epsilon c+\epsilon^2c_2+\cdots,\qquad
v_T=\epsilon d+\epsilon^2d_2+\cdots,
\]

\[
v_F=\epsilon^2b+\cdots,\qquad
v_H=\epsilon^2a+\cdots.
\]

Define

\[
q_1=s_t'\Lambda+Ro_tg_1,\qquad
A_0=-2\chi_ts_t,\qquad
C_0={\gamma s_t^2\over Ro_t},
\]

and

\[
A_1=-2\chi_tq_1+2\gamma s_t\theta_1,\qquad
C_1=\gamma\left({2s_tq_1\over Ro_t}
-{s_t^2\Lambda\over Ro_t^2}\right).
\]

The first grade is

\[
c''=0,\qquad d''=0,\qquad A_0c+C_0d=0.
\]

At the next grade, the previously missing meridional and radial equations
appear together:

\[
\boxed{
\begin{aligned}
a'+2b&=0,\\
b''+A_0c_2+C_0d_2&=-A_1c-C_1d,\\
c_2''+2\chi_ts_tb&=0,\\
d_2''&=0.
\end{aligned}}
\]

This calculation shows explicitly where the crossover parameter enters:
$q_1$ and $\theta_1$ are affine in $\Lambda$, so the forcing of the
second-grade radial equation is affine in $\Lambda$. No fold datum was used
to insert that dependence.

Let $(L,P,S,Q)$ multiply the four second-grade equations in the order shown.
The homogeneous second-grade formal adjoint is

\[
\begin{aligned}
-L'&=0,\\
2L+P''+2\chi_ts_tS&=0,\\
A_0P+S''&=0,\\
C_0P+Q''&=0,
\end{aligned}
\]

with inner concomitant

\[
\boxed{
\mathscr B_i=La+Pb'-P'b+Sc_2'-S'c_2+Qd_2'-Q'd_2.}
\]

The radial forcing contributes the inner candidate pairing

\[
-\int P(A_1c+C_1d)\,d\eta,
\]

which naturally separates into a constant part and a coefficient of
$\Lambda$. This is only the inner portion of $(\mathcal S_0,\mathcal S_1)$.
The outer correction and the overlap concomitant must still be included.

As a formal hierarchy bookkeeping device, let $(U,V,R,L,P,S,Q)$ multiply,
respectively,

\[
c'',\quad d'',\quad A_0c+C_0d,\quad
a'+2b,\quad b''+A_0c_2+C_0d_2+A_1c+C_1d,\quad
c_2''+2\chi_ts_tb,\quad d_2''.
\]

The full graded inner adjoint is then

\[
\begin{aligned}
U''+A_0R+A_1P&=0,\\
V''+C_0R+C_1P&=0,\\
-L'&=0,\\
2L+P''+2\chi_ts_tS&=0,\\
A_0P+S''&=0,\\
C_0P+Q''&=0,
\end{aligned}
\]

and the complete inner concomitant is

\[
\boxed{
\begin{aligned}
\mathscr B_i^{\rm graded}={}&Uc'-U'c+Vd'-V'd+La\\
&+Pb'-P'b+Sc_2'-S'c_2+Qd_2'-Q'd_2.
\end{aligned}}
\]

This block form records every integration-by-parts term, but the seven
multipliers are not seven independent physical adjoint fields. They must be
identified as successive coefficients of one four-component continuous
adjoint. The direct continuous derivation below supplies that identification
and supersedes any interpretation of this bookkeeping block as a standalone
adjoint BVP.

Matching the null mode to a regular outer Taylor expansion gives the
schematic large-$\eta$ conditions

\[
c\sim c_o'(0)\eta,\qquad d\sim d_o'(0)\eta,\qquad
a\sim a_o'(0)\eta,\qquad b\to b_o(0),
\]

with quadratic common parts in $c_2,d_2$. These polynomial common parts
must be subtracted consistently when comparing $\mathscr B_i$ with
$\mathscr B_o$; otherwise the apparent solvability integral depends on the
arbitrary overlap cutoff.

## Outer first correction of the null problem

Define the normalized continuous outer DAE for $X=(h,f,g,\theta)$ by

\[
\mathcal N(X;r)=
\begin{pmatrix}
h'+2f\\
\chi(s+rg)^2-\omega^2\\
g''+\chi\{r(2fg+hg')+2sf\}\\
\theta''+\mathrm{Pr}\,rh\theta'
\end{pmatrix}.
\]

At $(X_t,Ro_t)$, $\mathcal N_X=\mathcal L_o^{(0)}$. Separate the first
base-flow correction as

\[
X_1=X_1^{(0)}+\Lambda X_1^{(\mu)}.
\]

Here $X_1^{(0)}$ is generated by critical second-order inner matching data,
whereas $X_1^{(\mu)}$ contains the derivative of the generic outer family
with respect to $Ro$. Expanding the Jacobian gives

\[
\boxed{
\begin{aligned}
\mathcal L_o^{(0)}q_1=-{}&\mathcal N_{XX}[X_1^{(0)},q_0]\\
&-\Lambda\left(
\mathcal N_{XX}[X_1^{(\mu)},q_0]
+\mathcal N_{Xr}q_0\right).
\end{aligned}}
\]

Consequently,

\[
\mathcal L_o^{(1,0)}q_0=\mathcal N_{XX}[X_1^{(0)},q_0],
\]

\[
\mathcal L_o^{(1,\mu)}q_0=
\mathcal N_{XX}[X_1^{(\mu)},q_0]+\mathcal N_{Xr}q_0.
\]

There is no $O(\epsilon)$ contribution from the radial viscous outer term,
because that term is multiplied by $\epsilon^4$. Nullvector normalization
only fixes the multiple of $q_0$ contained in $q_1$ and cannot change the
Fredholm ratio.

## Primal common parts required by the concomitants

Write the leading outer null mode near $Y=0$ as

\[
\begin{aligned}
a_o(Y)&=A_1Y+\tfrac12A_2Y^2+\cdots,\\
b_o(Y)&=-\tfrac12A_1-\tfrac12A_2Y+\cdots,\\
c_o(Y)&=C_1Y+\tfrac12C_2Y^2+\cdots,\\
d_o(Y)&=D_1Y+\tfrac12D_2Y^2+\cdots.
\end{aligned}
\]

Continuity gives the expansion for $b_o$. The linearized algebraic
constraint gives

\[
2\chi_ts_tRo_tC_1-\gamma s_t^2D_1=0.
\]

Let $(a_{o1},b_{o1},c_{o1},d_{o1})$ be the components of $q_1$. Expanding
at $Y=\epsilon\eta$ and comparing with the inner hierarchy gives

\[
\begin{aligned}
c(\eta)&\sim C_1\eta,&
d(\eta)&\sim D_1\eta,\\
a(\eta)&\sim A_1\eta+a_{o1}(0),&
b(\eta)&\to-\tfrac12A_1,\\
c_2(\eta)&\sim\tfrac12C_2\eta^2+c_{o1}'(0)\eta+C_{20},&
d_2(\eta)&\sim\tfrac12D_2\eta^2+d_{o1}'(0)\eta+D_{20}.
\end{aligned}
\]

Matching also forces $c_{o1}(0)=d_{o1}(0)=0$. Otherwise an unmatched
constant would appear in the first-grade inner solutions $c''=d''=0$ with
homogeneous wall data. The constants $C_{20},D_{20}$ depend on the
higher-order expansion convention.

These are the common primal data to insert into
$\mathscr B_i^{\rm graded}$ and $\mathscr B_o$. Adjoint matching conditions
must follow by cancelling the coefficients of the independent common data.
Adding conditions directly at $Y=0$ would over-close the outer DAE.

## Exact four-field continuous adjoint

The graded calculation must be anchored to the adjoint of the original
four-field differential problem.

Write the physical residual as

\[
\begin{aligned}
R_H&=H'+2F,\\
R_F&=F''+\chi A_r+{\omega^2\over Ro},\\
R_G&=G''+\chi A_\theta,\\
R_T&=T''+\mathrm{Pr}\,Ro\,HT',
\end{aligned}
\]

where

\[
\chi=1-\gamma(T-1),\qquad q=s+Ro\,G,
\]

\[
A_r=Ro(F^2+HF')-{q^2\over Ro},\qquad
A_\theta=Ro(2FG+HG')+\mathrm{Co}\,F.
\]

For $v=(v_H,v_F,v_G,v_T)$, direct differentiation gives

\[
\begin{aligned}
(\mathcal Lv)_H={}&v_H'+2v_F,\\
(\mathcal Lv)_F={}&v_F''
+\chi RoF'v_H+\left(2\chi RoF+\chi RoH\partial_\eta\right)v_F
-2\chi qv_G-\gamma A_rv_T,\\
(\mathcal Lv)_G={}&v_G''
+\chi RoG'v_H+\chi(2RoG+\mathrm{Co})v_F\\
&+\left(2\chi RoF+\chi RoH\partial_\eta\right)v_G
-\gamma A_\theta v_T,\\
(\mathcal Lv)_T={}&v_T''
+\mathrm{Pr}\,Ro\,T'v_H+\mathrm{Pr}\,Ro\,H\partial_\eta v_T.
\end{aligned}
\]

Use the unweighted physical pairing
\(\int(\alpha R_H+\beta R_F+\zeta R_G+\tau R_T)\,d\eta\).  Term-by-term
integration by parts gives

\[
\int W\cdot\mathcal Lv\,d\eta
=\int (\mathcal L^\dagger W)\cdot v\,d\eta+
\left[\mathscr B_{\rm phys}(W,v)\right],
\qquad W=(\alpha,\beta,\zeta,\tau),
\]

where

\[
\begin{aligned}
(\mathcal L^\dagger W)_H={}&-\alpha'
+\chi RoF'\beta+\chi RoG'\zeta+\mathrm{Pr}\,RoT'\tau,\\
(\mathcal L^\dagger W)_F={}&2\alpha+\beta''
-(\chi RoH\beta)'+\chi(2RoG+\mathrm{Co})\zeta,\\
(\mathcal L^\dagger W)_G={}&-2\chi q\beta+\zeta''
-(\chi RoH\zeta)'+2\chi RoF\zeta,\\
(\mathcal L^\dagger W)_T={}&-\gamma A_r\beta-\gamma A_\theta\zeta+\tau''
-(\mathrm{Pr}\,RoH\tau)' ,
\end{aligned}
\]

and

\[
\boxed{
\begin{aligned}
\mathscr B_{\rm phys}={}&\alpha v_H
+\beta v_F'-\beta'v_F+\chi RoH\beta v_F\\
&+\zeta v_G'-\zeta'v_G+\chi RoH\zeta v_G\\
&+\tau v_T'-\tau'v_T+\mathrm{Pr}\,RoH\tau v_T.
\end{aligned}}
\]

There is no discretization or fold condition in these equations.  They are
the derivative of the basic-flow residual at a fixed member of the
crossover family.  A Euclidean left nullvector of a collocation Jacobian is
not automatically a discretization of $W$: its weights and normalization
must reproduce this pairing and boundary form.

### Regional power counting

In the outer variables take

\[
\alpha=\ell+O(\epsilon),\qquad
\beta=-\epsilon^2Ro_t\rho+O(\epsilon^3),\qquad
\zeta=\sigma+O(\epsilon),\qquad
\tau=\vartheta+O(\epsilon).
\]

The $\beta=O(\epsilon^2)$ scale is forced by the exact $G$-adjoint
equation: the algebraic term $-2\chi q\beta$ must balance
\(\zeta_{\eta\eta}=\epsilon^2\sigma_{YY}\).  The same scaling balances the
$\beta$- and $\vartheta_{\eta\eta}$-terms in the temperature adjoint.  The
$F$-adjoint balances $2\ell$ with the $O(1)$ azimuthal-adjoint term, so
\(\ell,\sigma,\vartheta=O(1)\).  Thus

\[
\boxed{
W_{\rm out}=(\alpha,\beta,\zeta,\tau)
=O(1,\epsilon^2,1,1).}
\]

The factor and sign in the physical trace are essential.  The outer DAE
uses the normalized radial residual \(\mathscr A=\chi q^2-\omega^2\),
whereas the leading physical radial residual is
\(R_F=-\mathscr A/Ro_t+O(\epsilon^4)\).  Equality of the two pairings
therefore gives

\[
\boxed{\beta_{\rm phys}=-\epsilon^2Ro_t\rho+O(\epsilon^3),}
\]

not \(\epsilon^2\rho\).  This conversion is required before comparing
physical radial traces with the DAE adjoint.

Substitution into the exact boundary form gives

\[
\mathscr B_{\rm phys}^{\rm out}
=\epsilon\mathscr B_o+O(\epsilon^2),
\]

where $\mathscr B_o$ is precisely the outer concomitant written above
(with $\vartheta$ in place of its earlier label $\tau$).  Within the
leading outer fields, the radial channel contributes only

\[
\left.\mathscr B_{\rm phys}^{\rm out}\right|_{\rm radial}
=-\epsilon^5Ro_t\{\rho b_Y-\rho_Yb+\chi Ro_t h\rho b\}.
\]

It therefore appears weak in the leading outer Green identity even though
it remains present in the exact four-field problem.

At fixed inner coordinate, use the primal scales already verified from the
fold mode.  The physical adjoint has the non-singular leading scale

\[
\alpha=L+O(\epsilon),\qquad
\beta=P+O(\epsilon),\qquad
\zeta=S+O(\epsilon),\qquad
\tau=Q+O(\epsilon),
\qquad
\boxed{W_{\rm in}=O(1,1,1,1).}
\]

Substitution into the exact adjoint equations at the first nontrivial inner
grade gives

\[
\begin{aligned}
-L'&=0,\\
2L+P''+2\chi_ts_tS&=0,\\
-2\chi_ts_tP+S''&=0,\\
{\gamma s_t^2\over Ro_t}P+Q''&=0.
\end{aligned}
\]

This is exactly the four-multiplier second-grade block above.  It identifies
$L,P,S,Q$ as the leading coefficients of the physical continuous adjoint.
The bookkeeping multipliers $(U,V,R)$ encode correction/constraint traces
introduced by retaining the first-grade equations after row weighting; they
are not three additional physical adjoint fields.

The inner limit of the exact boundary form begins as

\[
\begin{aligned}
\mathscr B_{\rm phys}^{\rm in}
={}&\epsilon\{Sc'-S'c+Qd'-Q'd\}\\
&+\epsilon^2\{La+Pb'-P'b+\text{first-correction terms from
\(\zeta,\tau\)}\}+O(\epsilon^3).
\end{aligned}
\]

The first-correction terms are represented by
$Uc'-U'c+Vd'-V'd$ in \(\mathscr B_i^{\rm graded}\).  Thus the
seven-multiplier concomitant is valid order-by-order bookkeeping, but only
the four-field expression defines the physical adjoint and normalization.

### Rule for deriving the adjoint interface conditions

Choose a movable overlap cut

\[
1\ll\eta_m\ll\epsilon^{-1},\qquad Y_m=\epsilon\eta_m.
\]

Applying the exact Green identity on $[0,\eta_m]$ and on
$[Y_m,\infty)$ shows that the interface term is

\[
\boxed{
\lim_{\substack{\eta_m\to\infty\\Y_m=\epsilon\eta_m\to0}}
\left[
\mathscr B_{\rm phys}^{\rm in}(\eta_m)
-\mathscr B_{\rm phys}^{\rm out}(Y_m)
\right]=0.}
\]

This is the source of the adjoint matching conditions.  It is not an
additional boundary condition imposed on the outer DAE.  Substitute the
primal common parts already listed above and expand both boundary forms in
the independent coefficients \(A_1,C_1,D_1,A_2,C_2,D_2,\ldots\).
Cancellation coefficient by coefficient determines the permitted entrance
traces of $(\ell,\rho,\sigma,\vartheta)$ and their relation to the
large-$\eta$ traces of $(L,P,S,Q)$ and their corrections.

Three consequences are immediate.

1. The leading inner boundary form contains the linearly growing $c,d$
   common part, and the next order contains the quadratic $c_2,d_2$ common
   part.  Both must be removed by common-part subtraction before a
   solvability integral is meaningful.
2. The outer entrance $Y=0$ is a matching boundary.  Conditions such as
   \(\sigma(0)=0\) or \(\vartheta(0)=0\) cannot be added merely to force the
   outer concomitant to vanish.
3. The Fredholm coefficients must be computed as a cutoff-independent sum,

\[
\mathcal S_j=
\mathcal S_j^{\rm in}(\eta_m)+\mathcal S_j^{\rm out}(Y_m),
\qquad j=0,1,
\]

with \(d\mathcal S_j/d\eta_m\to0\) in the overlap.  Only then is
\(\Lambda_f=-\mathcal S_0/\mathcal S_1\) a defined continuous selection
result.

The next calculation is now specific: solve the four-field adjoint in both
regions, use coefficient cancellation of the common parts to supply the
interface relations, and then test boundary-form cancellation, cutoff
invariance, adjoint-normalization invariance, and grid/mapping convergence
before quoting any numerical value of $\Lambda_f$.

## Canonical trace form at a movable interface

For bookkeeping, rewrite the exact physical concomitant as

\[
\mathscr B_{\rm phys}
=\alpha v_H+\beta v_F'
+\Pi_F v_F+\zeta v_G'+\Pi_G v_G
+\tau v_T'+\Pi_T v_T,
\]

where

\[
\Pi_F=-\beta'+\chi RoH\beta,\qquad
\Pi_G=-\zeta'+\chi RoH\zeta,\qquad
\Pi_T=-\tau'+\mathrm{Pr}\,RoH\tau.
\]

At a finite-$\epsilon$ artificial interface, the independent primal
canonical traces are
\(v_H,v_F,v_F',v_G,v_G',v_T,v_T'\).  Cancellation of the interface
concomitant therefore requires equality of their physical adjoint
coefficients:

\[
\boxed{
\alpha_i=\alpha_o,\quad
\beta_i=\beta_o,\quad
\Pi_{F,i}=\Pi_{F,o},\quad
\zeta_i=\zeta_o,\quad
\Pi_{G,i}=\Pi_{G,o},\quad
\tau_i=\tau_o,\quad
\Pi_{T,i}=\Pi_{T,o}.}
\]

These are not seven additional boundary conditions for the outer DAE.  They
are physical trace identities whose overlap expansions generate the regional
matching rules.  Since \(d/d\eta=\epsilon d/dY\) in the outer region and
\((\beta,v_F)=(O(\epsilon^2),O(\epsilon^2))\), the outer radial traces
enter at higher order than the $G,T$ traces.  One must therefore expand
these identities in the overlap and collect equal powers of $\epsilon$;
equating the leading inner coefficient $P$ directly to the outer coefficient
$\rho$ would be incorrect.

This mismatch means that the inner leading radial adjoint cannot be
identified directly with the outer radial multiplier.  A radial-adjoint
transition must be resolved by the overlap expansion before those
coefficients can be related.  A numerical implementation should monitor the
seven physical trace residuals above, not an unweighted difference between
rescaled inner and outer vectors.

## First overlap consequences and the radial-trace obstruction

Assume first that the outer adjoint has regular entrance traces and that the
leading inner azimuthal and thermal adjoint components remain bounded as
\(\eta\to\infty\).  The physical trace identities then immediately give

\[
\begin{aligned}
L_\infty&=\ell(0),&
S_\infty&=\sigma(0),&
Q_\infty&=\vartheta(0),\\
S'_\infty&=0,&
Q'_\infty&=0,&
P_\infty&=0.
\end{aligned}
\]

Here a subscript \(\infty\) denotes the large-\(\eta\) limit of the
leading inner coefficient.  The remaining leading relation follows from
either region:

\[
\boxed{
L_\infty=-\chi_t s_t S_\infty
=-\chi_t s_t\sigma(0).}
\]

Indeed, the outer $b$-adjoint equation at $Y=0$ is
\(2\ell(0)+2\chi_t s_t\sigma(0)=0\).  In the inner region, the same result
follows by taking the large-\(\eta\) limit of

\[
2L+P''+2\chi_t s_tS=0.
\]

Thus the $H$, $G$, and $T$ trace matching is already nontrivial but
consistent at leading order.  It is a direct consequence of the physical
concomitant, not a component-name identification.

The radial trace contains the unresolved singular part.  Let
\(R=S-S_\infty\).  The large-\(\eta\) leading inner equations reduce to

\[
P''=-2\chi_t s_tR,\qquad R''=2\chi_t s_tP,
\]

and hence

\[
\boxed{P''''+4\chi_t^2s_t^2P=0.}
\]

For the present critical state, \(\chi_ts_t>0\), so every bounded radial
inner mode is exponentially damped,

\[
P(\eta)=O\!\left(
e^{-\sqrt{\chi_ts_t}\,\eta}\right)
\quad\text{up to oscillatory factors}.
\]

By contrast, the physical outer radial trace is

\[
\beta_o=-\epsilon^2Ro_t\rho(0)+O(\epsilon^2Y,\epsilon^3),
\]

and its conjugate derivative trace is

\[
\Pi_{F,o}=
\epsilon^3\{Ro_t\rho_Y(0)-\chi Ro_t^2h(0)\rho(0)\}+\cdots.
\]

Consequently, the leading inner coefficient $P$ alone cannot supply the
outer $O(\epsilon^2)$ radial trace on an algebraic overlap
\(\eta=\epsilon^{-\alpha}\), \(0<\alpha<1\): its exponentially small tail
is beyond all algebraic orders there.  This does **not** yet imply
\(\rho(0)=0\).  It gives a precise dichotomy that must be resolved before
forming the Fredholm integrals:

1. a higher-grade inner adjoint coefficient produces the required
   $O(\epsilon^2)$ physical radial trace; or
2. an intermediate radial-adjoint transition layer connects the damped
   inner trace to the outer \(\rho(0)\) trace.

The exact physical expansion identifies the first candidate coefficient.
Write

\[
\begin{aligned}
\alpha&=L+\epsilon L_1+\cdots,\\
\beta&=P+\epsilon P_1+\epsilon^2P_2+\cdots,\\
\zeta&=S+\epsilon S_1+\cdots,\\
\tau&=Q+\epsilon Q_1+\cdots.
\end{aligned}
\]

Substitution in the four continuous adjoint equations gives the first
correction system

\[
\begin{aligned}
-L_1'+\chi_tRo_tg_1'S+\mathrm{Pr}\,Ro_t\theta_1'Q&=0,\\
2L_1+P_1''+2\chi_ts_tS_1-A_1S&=0,\\
A_0P_1+S_1''+A_1P&=0,\\
C_0P_1+Q_1''+C_1P&=0.
\end{aligned}
\]

This is obtained directly from the original four-field adjoint; in
particular, the same $A_1,C_1$ which carry $\Lambda$ in the graded primal
problem also enter the physical-adjoint correction.  It is therefore a
useful derivation check, not an added model assumption.

At the outer entrance,

\[
\beta_o=-\epsilon^2Ro_t\rho(0)+O(\epsilon^2Y,\epsilon^3).
\]

Consequently, a regular two-region matching expansion requires

\[
P\longrightarrow0,\qquad P_1\longrightarrow0,\qquad
P_2\longrightarrow-Ro_t\rho(0)
\quad\text{in the common-part sense}.
\]

Likewise \(P_2'\to0\) at this order, because the outer conjugate trace
\(\Pi_{F,o}\) starts at \(O(\epsilon^3)\).  Thus the first possible regular
inner provider of the outer radial trace is $P_2$, not $P$ or $P_1$.  The
next derivation is now sharper: derive the $P_2$ equation and its
large-\(\eta\) common part.  If that equation cannot admit the stated
constant trace, only then is an additional intermediate radial-adjoint
layer forced.  Until this test is completed, the relations for $L,S,Q$ are
valid leading interface data, whereas an entrance condition for $\rho$
remains deliberately unresolved.

## Second radial-trace grade: the local compatibility test

The preceding comparison must retain the normalization of the outer radial
row.  To test whether a regular two-region expansion is locally possible,
write, one order further,

\[
\begin{gathered}
Ro=Ro_t+\epsilon\Lambda+\epsilon^2r_2+\cdots,\qquad
G=\epsilon g_1+\epsilon^2g_2+\cdots,\\
T=T_t+\epsilon\theta_1+\epsilon^2\theta_2+\cdots,\qquad
F=\epsilon^2f_2+\cdots,\quad H=\epsilon^2h_2+\cdots,\\
\chi=\chi_t+\epsilon\chi_1+\epsilon^2\chi_2+\cdots,
\qquad \chi_1=-\gamma\theta_1,\quad \chi_2=-\gamma\theta_2.
\end{gathered}
\]

With \(s=s_t+\epsilon s_1+\epsilon^2s_2+\cdots\), define

\[
\begin{aligned}
q_2&=s_2+Ro_tg_2+\Lambda g_1,\\
A_2&=-2(\chi_tq_2+\chi_1q_1+\chi_2s_t),\\
C_2&=\gamma\left\{
 {q_1^2+2s_tq_2\over Ro_t}
 -{2s_tq_1\Lambda\over Ro_t^2}
 +s_t^2\left({\Lambda^2\over Ro_t^3}-{r_2\over Ro_t^2}\right)\right\},\\
\kappa_2&=\chi_tRo_t h_2.
\end{aligned}
\]

Here \(A_2\) is the $O(\epsilon^2)$ coefficient of \(-2\chi q\),
and \(-A_2\) is the corresponding coefficient of \(2\chi q\).  Direct
substitution into the exact four-field physical adjoint gives at
$O(\epsilon^2)$

\[
\begin{aligned}
-L_2'&+\chi_tRo_tf_2'P+\chi_tRo_tg_1'S_1
+\{(\chi_t\Lambda+\chi_1Ro_t)g_1'+\chi_tRo_tg_2'\}S\\
&\hspace{18mm}+\mathrm{Pr}\{Ro_t\theta_1'Q_1+
 (Ro_t\theta_2'+\Lambda\theta_1')Q\}=0,\\
2L_2&+P_2''-(\kappa_2P)'
 +2\chi_ts_tS_2-A_1S_1-A_2S=0,\\
A_0P_2&+A_1P_1+A_2P+S_2''-(\kappa_2S)'
 +2\chi_tRo_tf_2S=0,\\
C_0P_2&+C_1P_1+C_2P+Q_2''-2\gamma s_tf_2S
 -(\mathrm{Pr}\,Ro_th_2Q)'=0.
\end{aligned}
\]

This is the first grade at which the physical radial trace can be supplied
by a non-decaying inner coefficient.  It is not legitimate to discard the
terms containing \(h_2\) and \(f_2\): their common parts are the terms that
must reproduce the outer entrance relation.

For the critical base state, matching the base flow gives

\[
h_2(\eta)\sim h_1\eta+h_{20},\qquad
f_2(\eta)\to f_0,\qquad h_1+2f_0=0.
\]

Assume only the regular common-part conditions already required by the
physical traces, namely that \(P\) and \(P_1\) have no algebraic common
part and that \(S\to S_\infty\).  The third equation above then yields

\[
S_2''+A_0P_2\ \longrightarrow\
\chi_tRo_t(h_1-2f_0)S_\infty
=-4\chi_tRo_tf_0S_\infty.
\]

The corrected radial physical-trace identity gives

\[
\boxed{P_2\longrightarrow-Ro_t\rho(0).}
\]

Consequently,

\[
S_2''\longrightarrow
-2\chi_ts_tRo_t\rho(0)-4\chi_tRo_tf_0S_\infty.
\]

This is exactly the $Y\downarrow0$ limit of the outer $G$-adjoint equation:

\[
2\chi_ts_tRo_t\rho(0)+\sigma_{YY}(0)
-\chi_tRo_th_1\sigma(0)+2\chi_tRo_tf_0\sigma(0)=0,
\]

after using \(h_1=-2f_0\) and \(S_\infty=\sigma(0)\).  Hence a constant
second-grade radial trace is locally compatible with both regional
adjoint equations.  The earlier apparent sign mismatch was solely the
missing \(-Ro_t\) conversion between the normalized DAE multiplier
\(\rho\) and the physical multiplier \(\beta\).

The derivative trace first enters one grade later.  If
\(\beta=P+\epsilon P_1+\epsilon^2P_2+\epsilon^3P_3+\cdots\), then

\[
\Pi_{F,i}=-P'-\epsilon P_1'
+\epsilon^2\{-P_2'+\kappa_2P\}
+\epsilon^3\{-P_3'+\kappa_2P_1+\kappa_3P\}+\cdots.
\]

On the critical outer entrance, where \(h(0)=0\),

\[
\Pi_{F,o}=\epsilon^3Ro_t\rho_Y(0)+O(\epsilon^3Y,\epsilon^4).
\]

Thus regular common-part matching requires

\[
P_2'\to0,\qquad P_3'\to-Ro_t\rho_Y(0).
\]

The present calculation establishes only a local compatibility result:
the $P_2$ equation does not force \(\rho(0)=0\), so an intermediate
radial-adjoint layer is not forced by the first algebraic overlap test.
It does not prove that the global $P_2,P_3$ inner problem and the outer DAE
adjoint admit compatible boundary data.  That BVP, followed by direct
concomitant and cutoff-independence tests, remains necessary before forming
\(\mathcal S_0,\mathcal S_1\) or quoting \(\Lambda_f\).

## Seven physical traces: common-part ledger

The interface must be formulated in the canonical physical trace vector,
not in the four multiplier names used by either regional system.  Define

\[
\mathbf Z=(\alpha,\beta,\Pi_F,\zeta,\Pi_G,\tau,\Pi_T),
\qquad
\mathbf t=(v_H,v_F',v_F,v_G',v_G,v_T',v_T),
\]

so that \(\mathscr B_{\rm phys}=\mathbf Z\mathbin{\cdot}\mathbf t\).
At a movable cut, the required statement is equality of the *common parts*
of these seven coefficients.  It is not an equality of a raw inner vector
with a raw outer vector.

Write the regular outer entrance expansions as

\[
\begin{gathered}
\ell=\ell_0+Y\ell_{Y0}+\cdots,\quad
\rho=\rho_0+Y\rho_{Y0}+\cdots,\quad
\sigma=\sigma_0+Y\sigma_{Y0}+\tfrac12Y^2\sigma_{YY0}+\cdots,\\
\vartheta=\vartheta_0+Y\vartheta_{Y0}
 +\tfrac12Y^2\vartheta_{YY0}+\cdots,\qquad
h=h_1Y+\cdots .
\end{gathered}
\]

Here \(h(0)=0\) is the critical outer entrance condition, while the
coefficient \(h_1\) is the same coefficient that appears in the inner
common part \(h_2\sim h_1\eta+h_{20}\).  The seven outer physical traces,
through their first nonzero common parts, are

\[
\begin{aligned}
\alpha_o&=\ell_0+Y\ell_{Y0}+\cdots,\\
\beta_o&=-\epsilon^2Ro_t\{\rho_0+Y\rho_{Y0}+\cdots\}+O(\epsilon^3),\\
\Pi_{F,o}&=\epsilon^3Ro_t\rho_{Y0}+O(\epsilon^3Y,\epsilon^4),\\
\zeta_o&=\sigma_0+Y\sigma_{Y0}+\tfrac12Y^2\sigma_{YY0}+\cdots,\\
\Pi_{G,o}&=\epsilon\{-\sigma_{Y0}
 +Y(-\sigma_{YY0}+\chi_tRo_th_1\sigma_0)+\cdots\}+O(\epsilon^2),\\
\tau_o&=\vartheta_0+Y\vartheta_{Y0}
 +\tfrac12Y^2\vartheta_{YY0}+\cdots,\\
\Pi_{T,o}&=\epsilon\{-\vartheta_{Y0}
 +Y(-\vartheta_{YY0}+\mathrm{Pr}\,Ro_th_1\vartheta_0)+\cdots\}
 +O(\epsilon^2).
\end{aligned}
\]

The \(O(\epsilon^3Y)\) radial term contains
\(-\chi_tRo_t^2h\rho\).  It vanishes at the critical entrance but must be
retained when the next radial grade is derived.  This is why an entrance
comparison of \(\beta\) alone is insufficient.

On the inner side, write

\[
\begin{aligned}
\alpha_i&=L+\epsilon L_1+\cdots,\\
\beta_i&=P+\epsilon P_1+\epsilon^2P_2+\epsilon^3P_3+\cdots,\\
\zeta_i&=S+\epsilon S_1+\epsilon^2S_2+\cdots,\\
\tau_i&=Q+\epsilon Q_1+\epsilon^2Q_2+\cdots .
\end{aligned}
\]

Substituting \(Y=\epsilon\eta\), cancellation of the seven independent
primal trace coefficients gives the following common-part ledger:

\[
\begin{aligned}
L&\to\ell_0,&
L_1&\sim\eta\ell_{Y0}+\ell_{10},\\
P&\to0,&P_1&\to0,&P_2&\to-Ro_t\rho_0,&P_2'&\to0,&P_3'&\to-Ro_t\rho_{Y0},\\
S&\to\sigma_0,&S_1&\sim\eta\sigma_{Y0}+\sigma_{10},&
S_2&\sim\tfrac12\eta^2\sigma_{YY0}+\eta\sigma_{21}+\sigma_{20},\\
Q&\to\vartheta_0,&Q_1&\sim\eta\vartheta_{Y0}+\vartheta_{10},&
Q_2&\sim\tfrac12\eta^2\vartheta_{YY0}+\eta\vartheta_{21}+\vartheta_{20}.
\end{aligned}
\]

The constants with second subscripts are convention-dependent homogeneous
common parts; they are not independent physical interface data.  In
particular, this ledger recovers the already derived leading conditions
\(S'\to Q'\to0\), rather than imposing them separately.

The conjugate derivative traces provide a direct check that the ledger is
consistent.  Since \(h_2\sim h_1\eta+h_{20}\),

\[
\begin{aligned}
\Pi_{G,i}
 &=-S'-\epsilon S_1'
 +\epsilon^2\{-S_2'+\chi_tRo_th_2S\}+\cdots,\\
\Pi_{T,i}
 &=-Q'-\epsilon Q_1'
 +\epsilon^2\{-Q_2'+\mathrm{Pr}\,Ro_th_2Q\}+\cdots.
\end{aligned}
\]

Their common parts become, respectively,

\[
\begin{aligned}
\Pi_{G,i}&=-\epsilon\sigma_{Y0}
 +\epsilon^2\eta(-\sigma_{YY0}+\chi_tRo_th_1\sigma_0)+\cdots,\\
\Pi_{T,i}&=-\epsilon\vartheta_{Y0}
 +\epsilon^2\eta(-\vartheta_{YY0}
 +\mathrm{Pr}\,Ro_th_1\vartheta_0)+\cdots,
\end{aligned}
\]

which are exactly \(\Pi_{G,o}\) and \(\Pi_{T,o}\) after \(Y=\epsilon\eta\).
The $H$ trace is also compatible: the inner $L_1$ equation produces

\[
L_1'\to\chi_tRo_tg_Y(0)\sigma_0+
\mathrm{Pr}\,Ro_t\theta_Y(0)\vartheta_0=\ell_{Y0},
\]

which is the $Y=0$ outer $a$-adjoint equation.  Thus the first nontrivial
common parts of the $H$, $G$, and $T$ channels cancel without an imposed
outer boundary condition.

The two second-derivative entrance equations give a stronger local check.
The $G$ channel gives

\[
\sigma_{YY0}=-2\chi_ts_tRo_t\rho_0
 +\chi_tRo_th_1\sigma_0-2\chi_tRo_tf_0\sigma_0,
\]

which is the $P_2$ compatibility calculation above.  The inner thermal
adjoint equation similarly gives

\[
Q_2''\to\gamma s_t^2\rho_0+2\gamma s_tf_0\sigma_0
 +\mathrm{Pr}\,Ro_th_1\vartheta_0,
\]

identical to the outer $T$-adjoint entrance equation after using
\(B(0)=2s_tf_0\).  These two equalities verify the first algebraic
common-part cancellation of every channel for which an entrance derivative
is determined locally.

This is still not a complete adjoint BVP.  The radial values
\(\rho_0,\rho_{Y0}\), the homogeneous common-part constants, and the
far-field adjoint conditions must be selected by the coupled inner--outer
problem.  The next executable step is therefore to solve that problem with
the ledger above as its interface map, then verify

\[
\mathscr B_i(\eta_m)-\mathscr B_o(\epsilon\eta_m)\to0,
\qquad
{d\over d\eta_m}
\left[\mathcal S_j^{\rm in}(\eta_m)+
\mathcal S_j^{\rm out}(\epsilon\eta_m)\right]\to0
\quad(j=0,1),
\]

over a range of overlap cuts.  Only those cutoff-independent sums define
the continuous Fredholm coefficients.

## Third radial-trace grade

The value trace of the radial multiplier is supplied by \(P_2\), but its
conjugate derivative trace first appears one grade later.  To make that
statement operational, define the coefficients of the exact products

\[
\chi q=K_0+\epsilon K_1+\epsilon^2K_2+\epsilon^3K_3+\cdots,
\qquad
\chi RoH=\epsilon^2\kappa_2+\epsilon^3\kappa_3+\cdots,
\]

where

\[
\begin{aligned}
K_0&=\chi_ts_t,\\
K_1&=\chi_tq_1+\chi_1s_t,\\
K_2&=\chi_tq_2+\chi_1q_1+\chi_2s_t,\\
K_3&=\chi_tq_3+\chi_1q_2+\chi_2q_1+\chi_3s_t,\\
\kappa_3&=\chi_tRo_th_3+(\chi_t\Lambda+\chi_1Ro_t)h_2.
\end{aligned}
\]

Similarly write

\[
\chi RoF=\epsilon^2\mu_2+\epsilon^3\mu_3+\cdots,
\qquad
\mu_2=\chi_tRo_tf_2,
\quad
\mu_3=\chi_tRo_tf_3+(\chi_t\Lambda+\chi_1Ro_t)f_2.
\]

The \(O(\epsilon^3)\) coefficient of the exact \(G\)-adjoint equation is

\[
\boxed{
\begin{aligned}
0={}&S_3''-2K_0P_3-2K_1P_2-2K_2P_1-2K_3P\\
&-(\kappa_2S_1+\kappa_3S)'
 +2\mu_2S_1+2\mu_3S .
\end{aligned}}
\]

The corresponding \(F\)-adjoint equation contains the radial second
derivative at the same grade:

\[
\boxed{
2L_3+P_3''-(\kappa_2P_1+\kappa_3P)'
 +2\chi_ts_tS_3-A_1S_2-A_2S_1-A_3S=0.}
\]

where \(A_3\) is the \(O(\epsilon^3)\) coefficient defined by
\(2\chi q=2\chi_ts_t-A_1\epsilon-A_2\epsilon^2-A_3\epsilon^3+\cdots\).
These equations show exactly where the first radial derivative matching data
enter; no new physical multiplier is introduced.

The physical radial conjugate trace is

\[
\Pi_{F,i}=-P'-\epsilon P_1'
 +\epsilon^2(-P_2'+\kappa_2P)
 +\epsilon^3(-P_3'+\kappa_2P_1+\kappa_3P)+\cdots,
\]

whereas the outer trace is

\[
\Pi_{F,o}=\epsilon^3Ro_t\rho_{Y0}+O(\epsilon^3Y,\epsilon^4)
\]

at the critical entrance.  Since \(P,P_1\) have no algebraic common part,
the common-part conditions are therefore

\[
P_2'\to0,
\qquad
P_3'\to-Ro_t\rho_{Y0}.
\]

The second condition means that \(P_3\) generally has a linear common part,
\(P_3\sim-Ro_t\rho_{Y0}\eta+\text{constant}\), which is precisely the
inner representation of the outer expansion
\(-\epsilon^2Ro_t\rho(\epsilon\eta)\).  It does not impose
\(\rho_{Y0}=0\).  A nonzero \(\rho_{Y0}\) is compatible with the overlap
provided the \(S_3\) and \(P_3\) equations reproduce the same polynomial
common part.  Checking that reproduction is the next local algebraic test;
global existence and cutoff independence still require solving the coupled
regional BVP.

## Explicit \(Y\)-derivative entrance test

The first test involving \(\rho_{Y0}\) follows directly from the leading
outer \(G\)-adjoint equation.  Define \(K=\chi q\), \(M=\chi Ro_th\), and
\(N=\chi Ro_tf\), so that

\[
2Ro_tK\rho+\sigma_{YY}-(M\sigma)_Y+2N\sigma=0.
\]

At the critical entrance \(M(0)=0\).  Differentiating once in \(Y\) gives

\[
\boxed{
\sigma_{YYY0}
 +2Ro_t(K_{Y0}\rho_0+K_0\rho_{Y0})
 -(M_{YY0}\sigma_0+2M_{Y0}\sigma_{Y0})
 +2(N_{Y0}\sigma_0+N_0\sigma_{Y0})=0.}
\]

This is the outer target for the linear common part of the third-grade
inner equation: \(S_3''\sim\sigma_{YYY0}\eta+\cdots\) and
\(P_3\sim-Ro_t\rho_{Y0}\eta+\cdots\).  After all remaining third-grade
base and adjoint coefficients are expanded, the coefficient of \(\eta\) in
the \(S_3\) equation must reduce to the boxed identity.  The test checks
local compatibility of a nonzero \(\rho_{Y0}\); it does not determine that
quantity or replace the global matched-adjoint BVP.

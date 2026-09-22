# Generic thermal-tail crossover study (2026-09-09)

## Question and separation from the fold calculation

Let

\[
\mu=Ro-Ro_t,
\qquad Ro_t=-0.470367416464038.
\]

The generic outer/inner calculations show, locally,

\[
h_0\sim c_h\mu,
\qquad D_1\sim c_D\mu,
\qquad \|(F_1,H_1)\|\sim c_M|\mu|,
\]

with $c_h=-7.785373304667$, $c_D=1.719419134040$, and
$c_M=19.06144962830$. Therefore the generic contribution
$\epsilon(F_1,H_1)$ becomes the same order as the critical contribution
$\epsilon^2(F_2,H_2)$ when

\[
\mu=O(\epsilon).
\]

The distinguished variable to test is consequently

\[
\Lambda={Ro-Ro_t\over\epsilon}.
\]

This study uses only the basic-flow boundary-value problem with prescribed
$(Ro,H_\infty)$, releases $T_w$, and imposes no fold or null-vector
condition. Corrected finite-$\epsilon$ fold data will be compared only
after this no-fold test has been completed.

## Planned no-fold sweep

For each fixed

\[
\Lambda\in\{-0.4,-0.3,0,0.5,1,2\},
\qquad
\epsilon\in\{0.03,0.02,0.015,0.01\},
\]

solve at

\[
Ro=Ro_t+\epsilon\Lambda,
\qquad H_\infty=-\epsilon.
\]

Use mapped collocation degrees $N=240$ and $N=320$, with the established
mapping $(a,b,c)=(2,0.6,0.5)$, $Pr=0.72$, and $\gamma=1$.
The regenerated data are written under
`work/results/generic_tail_crossover_20260909/`; preserved baseline data are
not modified.

The primary falsifiable prediction is fixed-$\Lambda$ collapse of

\[
{F\over\epsilon^2},\qquad {H\over\epsilon^2},\qquad
{G\over\epsilon},\qquad {T-T_t\over\epsilon}
\]

on every fixed inner interval.  In addition,

\[
{h_0(Ro)\over\epsilon}\longrightarrow c_h\Lambda
\]

provides an independently calculated outer-entry check.  Failure of these
collapses under grid refinement would reject, or restrict, the proposed
distinguished crossover scaling.

## Local two-parameter derivation

Set

\[
Ro=Ro_t+\epsilon\Lambda,
\]

where $\Lambda=O(1)$. The generic fixed-$Ro$ expansion has

\[
F=\epsilon F_1(\eta;Ro)+\epsilon^2F_2+\cdots,
\qquad
H=\epsilon H_1(\eta;Ro)+\epsilon^2H_2+\cdots.
\]

The independently solved generic first-order inner problem gives
$F_1(\eta;Ro_t)=H_1(\eta;Ro_t)=0$ and the local outer-family calculation
shows differentiable, nonzero first derivatives with respect to $Ro$.
Taylor expansion in the second small parameter therefore gives

\[
\begin{aligned}
F_1(\eta;Ro_t+\epsilon\Lambda)
 &=\epsilon\Lambda\,\partial_{Ro}F_1(\eta;Ro_t)+O(\epsilon^2),\\
H_1(\eta;Ro_t+\epsilon\Lambda)
 &=\epsilon\Lambda\,\partial_{Ro}H_1(\eta;Ro_t)+O(\epsilon^2).
\end{aligned}
\]

Consequently the uniform local form is

\[
\boxed{
\begin{aligned}
F&=\epsilon^2\left[F_2^{(t)}
 +\Lambda\,\partial_{Ro}F_1|_t\right]+o(\epsilon^2),\\
H&=\epsilon^2\left[H_2^{(t)}
 +\Lambda\,\partial_{Ro}H_1|_t\right]+o(\epsilon^2).
\end{aligned}}
\]

Thus the critical $O(\epsilon^2)$ fields and the vanishing generic
$O(\epsilon)$ fields do not compete until $\mu/\epsilon=\Lambda=O(1)$.
This is the direct origin of the distinguished scaling; it is not inferred
from a fold trajectory.

For $G$ and $T$, write

\[
G=\epsilon g_1+\cdots,
\qquad
T=T_t+\epsilon\theta_1+\cdots.
\]

The leading inner equations give

\[
g_1=p_t\eta,
\qquad
\theta_1=\mathcal T'(Ro_t)\Lambda+B_T\eta.
\]

Substitution into the $O(\epsilon)$ radial algebraic balance yields

\[
\gamma s_t^2\theta_1
-2\chi_t s_tRo_tg_1+D_t'\Lambda=0.
\]

Its constant part vanishes because
$\mathcal T'(Ro_t)=-D_t'/(\gamma s_t^2)$, and its coefficient of
$\eta$ vanishes by the already established entrance-slope relation.
Hence no additional fold constraint appears at this order: the crossover
ansatz is internally compatible for arbitrary bounded $\Lambda$.

The leading outer problem is likewise still the critical outer solution.
The parameter $\Lambda$ first appears in its linear correction, whose
entrance data include

\[
h_1^{\mathrm{out}}(0)=c_h\Lambda.
\]

Accordingly, a complete analytic crossover BVP would couple the second-order
inner problem to the first correction of the critical outer problem.  The
fixed-$\Lambda$ numerical sweep tests the scaling before undertaking that
larger derivation.

## Numerical results

All 48 no-fold BVPs converged (six values of $\Lambda$, four values of
$\epsilon$, and two grids). The largest reported nonlinear residual was
$3.2\times10^{-9}$. On the sampled fixed inner interval
$\eta=0.5,1,2,5,8$, comparison with the $\epsilon=0.01$ scaled profile gives
the following worst-case relative differences over all tested $\Lambda$ on
the $N=320$ grid:

| field | $\epsilon=0.03$ | $\epsilon=0.02$ | $\epsilon=0.015$ |
|---|---:|---:|---:|
| $H/\epsilon^2$ | 0.192 | 0.103 | 0.0534 |
| $F/\epsilon^2$ | 0.130 | 0.0689 | 0.0355 |
| $G/\epsilon$ | 0.124 | 0.0643 | 0.0327 |
| $(T-T_t)/\epsilon$ | 0.0255 | 0.0127 | 0.00633 |

The differences decrease monotonically as $\epsilon$ decreases. At
$\epsilon=0.01$, the worst $N=240$ versus $N=320$ relative difference is
0.00211 for $F/\epsilon^2$, 0.000629 for $H/\epsilon^2$, and below
$2.3\times10^{-5}$ for the other two scaled fields. The collapse is thus
well separated from the grid error.

The algebraic wall defect supplies a second test. A fit of
$D/\epsilon^2=a(\epsilon)+b(\epsilon)\Lambda$ over the full sampled interval
gives, on $N=320$,

| $\epsilon$ | $a$ | $b$ | fit RMS |
|---:|---:|---:|---:|
| 0.030 | 0.257863 | 1.64809 | 0.0404 |
| 0.020 | 0.256674 | 1.68116 | 0.0234 |
| 0.015 | 0.255885 | 1.69444 | 0.0161 |
| 0.010 | 0.254729 | 1.70536 | 0.00972 |

The slope approaches the independently obtained generic coefficient
$c_D=1.719419$, while the nonzero intercept is the critical second-order
defect. Both the field collapse and this affine defect law support the
distinguished crossover.

## Fold comparison after the no-fold test

Using the already corrected fold chain only after establishing the basic-flow
crossover gives

| $\epsilon$ | $\Lambda_f$, $N=280$ | $\Lambda_f$, $N=320$ |
|---:|---:|---:|
| 0.10 | -0.339119 | -0.339119 |
| 0.05 | -0.316586 | -0.316586 |
| 0.02 | -0.302896 | -0.302897 |
| 0.01 | -0.295249 | -0.298372 |
| 0.005 | -0.365414 | -0.344199 |

The resolved range $0.01\leq\epsilon\leq0.1$ has bounded
$\Lambda_f=O(1)$ and trends toward approximately $-0.30$, so the finite-
$\epsilon$ fold does enter the distinguished crossover layer. The
$\epsilon=0.005$ values are not grid converged and are excluded from any
limit estimate. These data do not yet derive the selected limiting value
of $\Lambda_f$; that requires the next-order augmented or equivalent
fold-solvability problem.

## Evidence-level conclusion

The no-fold calculations independently establish that the generic
$O(\epsilon)$ meridional contribution and the critical $O(\epsilon^2)$
contribution merge under $Ro-Ro_t=\epsilon\Lambda$. Corrected fold data then
show that the fold trajectory enters this layer. What remains unproved is
the full coupled second-order-inner/first-correction-outer BVP and the
next-order solvability condition that selects the fold path within the
layer.

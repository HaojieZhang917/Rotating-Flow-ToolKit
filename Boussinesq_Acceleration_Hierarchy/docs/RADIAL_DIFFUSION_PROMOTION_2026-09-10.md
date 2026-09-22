# Radial-diffusion promotion at the heated-BEK fold--tail limit

## Scope

This note defines Gate R4. It restores the radial viscous and thermal
diffusion omitted from the parabolized non-similar model, projects the
resulting local operators onto the finite-epsilon fold null pair, and tests
their scaling. This is a frozen-radius, steady, axisymmetric diagnostic. It
is not a finite-disk solution, temporal stability calculation, or full
three-dimensional stability result.

## Cylindrical diffusion identities

For an axisymmetric scalar (q(R,z)),

\[
 \nabla^2q=q_{RR}+R^{-1}q_R+q_{zz}.
\]

For (u_R=RA(R,z)) and (u_\phi=RB(R,z)), the radial and azimuthal vector
Laplacian components are

\[
 (\nabla^2\boldsymbol u)_R
 =R A_{RR}+3A_R+RA_{zz},
\]

\[
 (\nabla^2\boldsymbol u)_\phi
 =R B_{RR}+3B_R+RB_{zz}.
\]

For the axial velocity (W), the viscous expression would be
(W_{RR}+R^{-1}W_R+W_{zz}). It is not inserted into the present (H)
row: that row is continuity, not axial momentum. Retaining axial viscous
diffusion requires adding the axial momentum equation and its pressure
correction. Adding an (H)-diffusion term to continuity would be
incorrect.

## Diffusive non-similar BEK system

With (x=R/L), (\eta=z/\ell), (\delta_R=\ell/L), and the definitions
in `NONSIMILAR_HEATED_BEK_BOUNDARY_LAYER_2026-09-09.md`, the three transport
residuals acquire

\[
 \delta_R^2(F_{xx}+3x^{-1}F_x),\qquad
 \delta_R^2(G_{xx}+3x^{-1}G_x),\qquad
 \delta_R^2(T_{xx}+x^{-1}T_x),
\]

respectively. Continuity is unchanged. Setting (\xi=\log x) gives

\[
 F_{xx}+3x^{-1}F_x=x^{-2}(F_{\xi\xi}+2F_\xi),
\]

\[
 G_{xx}+3x^{-1}G_x=x^{-2}(G_{\xi\xi}+2G_\xi),
\qquad
 T_{xx}+x^{-1}T_x=x^{-2}T_{\xi\xi}.
\]

At a frozen physical radius (R_0=Lx_0), define

\[
 \kappa_0={\delta_R^2\over x_0^2}=left({\ell\over R_0}\right)^2.
\]

The local Mellin pencil is therefore

\[
 \boxed{J+\lambda M_r+\kappa_0(\lambda^2D_2+\lambda D_1)}.
\]

In field order ((H,F,G,T)), the nonzero interior blocks are

\[
 (D_2)_{FF}=(D_2)_{GG}=(D_2)_{TT}=I,
\]

\[
 (D_1)_{FF}=(D_1)_{GG}=2I.
\]

All replacement boundary rows vanish. A Julia algebraic regression checks
these blocks for five real (\lambda) values with zero discrepancy.

## Fold projection and normalization

Use the right-mode convention (\max|v_G|=1), fixed sign, and divide all
Fredholm projections by the ordinary wall-temperature transversality:

\[
 \widehat c_{rr}={w^TD_2v\over w^TR_{T_w}},\qquad
 \widehat c_{r,{\rm geo}}={w^TD_1v\over w^TR_{T_w}}.
\]

This avoids the artificial poles produced when the discrete (w^Tv)
passes close to zero. It is also invariant under rescaling of (w).

## Numerical result

The corrected-energy fold profiles give:

| ϵ | ĉ_rr, N=280 | ĉ_rr, N=320 | ĉ_r,geo, N=280 | ĉ_r,geo, N=320 |
|---:|---:|---:|---:|---:|
| 0.10 | -48.4305 | -48.4321 | 16.8143 | 16.8149 |
| 0.05 | -97.5413 | -97.5448 | 31.2143 | 31.2153 |
| 0.02 | -244.687 | -243.619 | 74.7442 | 74.6311 |
| 0.01 | -505.996 | -479.476 | 148.687 | 146.570 |

The ϵ=0.005 point is retained in the raw table but excluded from the
fit because the saved fold and the (c_2) projection are not fully
resolved there. Four-point fits over ϵ=0.10--0.01 give

\[
 |\widehat c_{rr}|\simeq4.64\epsilon^{-1.017}\quad(N=280),
\]

\[
 |\widehat c_{rr}|\simeq4.91\epsilon^{-0.996}\quad(N=320),
\]

while the cylindrical first-derivative projection has powers (-0.947)
and (-0.942). Hence the proposed (O(1)) diffusion coefficient is
falsified:

\[
 \boxed{\widehat c_{rr}=O(\epsilon^{-1}),\qquad
 \widehat c_{r,{\rm geo}}=O(\epsilon^{-1}).}
\]

At (N=320), the decomposition of ĉ_rr into (F,G,T) blocks shows that
the divergence is dominated by the thermal component. For ϵ=0.10 and
0.01, respectively,

\[
 (\widehat c_{rr}^{F},\widehat c_{rr}^{G},\widehat c_{rr}^{T})
 =(-0.0135,8.421,-56.840),
\]

\[
 (-0.000125,73.285,-552.761).
\]

The (G) contribution partially cancels, but does not remove, the thermal
growth.

## Promoted radial scale

The previously verified first-order radial coefficient satisfies
(\widehat c_r=O(\epsilon)). For variations on an (O(1)) scale in
(\xi), the projected radial terms have orders

\[
 \widehat c_rA_\xi=O(\epsilon),
\]

\[
 \kappa_0(\widehat c_{rr}A_{\xi\xi}
 +\widehat c_{r,{\rm geo}}A_\xi)
 =O(\kappa_0\epsilon^{-1}).
\]

They balance when

\[
 \kappa_0=O(\epsilon^2),
 \qquad
 \boxed{{R_0\over\ell}=O(\epsilon^{-1}).}
\]

Thus the initially proposed (R_0/\ell=O(\epsilon^{-1/2})) scaling would
follow only if ĉ_rr were (O(1)); the computed adjoint projection shows
that it is not. The thermal-tail response supplies an additional
(\epsilon^{-1}) amplification and moves the candidate interaction radius
to (O(\ell/\epsilon)).

This balance assumes (A) varies on an (O(1)) logarithmic-radial scale.
A separate slow or rapid (\xi) scale changes the derivative powers and
must be derived in the subsequent multiple-scale problem. Furthermore,
(\kappa_0=e^{-2\xi}\delta_R^2) prevents a global constant-coefficient
Mellin interpretation.

## Evidence boundary and next Gate

An independent map check recomputed the (N=240) fold chain, rather than
only reprojecting old profiles, at (b=0.55) and (b=0.65), keeping
((a,c)=(2,0.5)). The four-point exponents are

| map (b) | power of ĉ_r | power of ĉ_rr | power of ĉ_r,geo |
|---:|---:|---:|---:|
| 0.55 | 1.03044 | -1.00086 | -0.94336 |
| 0.60, N=320 | 1.03079 | -0.99608 | -0.94184 |
| 0.65 | 1.03095 | -1.01281 | -0.94629 |

Thus the powers relevant to promotion are insensitive to this rational-map
variation. The prefactor of ĉ_rr remains approximately (4.7)--(4.9).
Gate R4 therefore passes both the grid and first independent mapping tests.

The next task is to derive a variable-radius amplitude equation retaining
both cylindrical diffusion terms; do not start a full finite-disk
calculation until that reduced problem determines the necessary radial
domain and edge data.

Sources and outputs:

- `work/src/BEKRadialResponse.jl`;
- `work/scripts/verify_radial_diffusion_operator_20260910.jl`;
- `work/scripts/analyze_fold_radial_transversality_20260909.jl`;
- `work/results/fold_radial_transversality_20260909/`.
- `work/results/corrected_energy_epsilon_fold_chain/map_b055_N240/`;
- `work/results/corrected_energy_epsilon_fold_chain/map_b065_N240/`.

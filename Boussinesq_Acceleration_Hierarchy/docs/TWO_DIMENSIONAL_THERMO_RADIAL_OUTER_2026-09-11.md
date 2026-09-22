# Distinguished two-dimensional thermo-radial outer region

## Scaling

Gate R4 gives

\[
\widehat c_r=O(\epsilon),\qquad
\widehat c_{rr}=O(\epsilon^{-1}),
\]

so radial diffusion meets radial transport at \(R/\ell=O(\epsilon^{-1})\).
The wall-normal thermal tail already has \(z/\ell=O(\epsilon^{-1})\).
Introduce

\[
X=\epsilon R/\ell,\qquad Y=\epsilon z/\ell=\epsilon\eta,
\]

and the critical outer fields

\[
F=\epsilon^2f(X,Y),\quad H=\epsilon h(X,Y),\quad
G=g(X,Y),\quad T=t(X,Y).
\]

## Leading thermo-azimuthal subsystem

Let \(q=s_t+Ro_tg\), \(\chi=1-\gamma(t-1)\),
\(Co_t=2s_t\), and \(\omega_t=s_t+Ro_t\). Direct substitution into the
radially diffusive four-field equations gives

\[
\boxed{h_Y+2f+Xf_X=0,}
\]

\[
\boxed{\chi q^2=\omega_t^2,}
\]

\[
\boxed{
g_{YY}+g_{XX}+\frac3Xg_X+
\chi\left[Ro_t(2fg+Xfg_X+hg_Y)+Co_tf\right]=0,
}
\]

\[
\boxed{
t_{YY}+t_{XX}+\frac1Xt_X+
Pr\,Ro_t(Xft_X+ht_Y)=0.
}
\]

Setting all \(X\)-derivatives to zero recovers the existing one-dimensional
critical outer DAE algebraically. Thus the promoted region is genuinely
two-dimensional in its thermal and azimuthal diffusion. The subsequent
full-axisymmetric audit shows that these equations are necessary but not a
complete full-NS closure.

The leading radial centrifugal balance remains algebraic. Meridional
inertia, \(f_{YY}\), radial \(f\)-diffusion, full axial momentum,
\(W\)-viscosity, and the pressure correction enter four powers later in
that balance. However, this is their first nonzero meridional momentum
order, and pressure compatibility supplies an additional constraint on
\(f,h\). In particular, \(H\)-diffusion must not be inserted into
continuity. See FULL_AXISYMMETRIC_OUTER_HIERARCHY_AUDIT_2026-09-11.md.

## Full-equation power count

At \(R,z=O(\ell/\epsilon)\),

| balance | physical order |
|---|---:|
| centrifugal acceleration and leading radial pressure | \(\Omega^2\ell/\epsilon\) |
| azimuthal advection and two-directional viscosity | \(\Omega^2\ell\epsilon\) |
| thermal advection and two-directional diffusion | \(\Omega\epsilon^2\) |
| meridional inertia and viscosity | \(\Omega^2\ell\epsilon^3\) |

Leading axial momentum makes the rotational pressure independent of \(Y\);
the uniform rigid far field then yields \(\chi q^2=\omega_t^2\). The first
nonzero axial/radial meridional subsystem determines the pressure
correction and imposes a compatibility condition on the displayed fields;
it is therefore part of the full distinguished hierarchy even though its
absolute momentum order is smaller.

## Why the diffusion projection grows

For \(O(1)\) outer direct and adjoint thermal modes,

\[
\int_0^\infty w_TD_2v_T\,d\eta
=\frac1\epsilon\int_0^\infty
w_T^{(o)}D_2v_T^{(o)}\,dY.
\]

The outer integration length therefore supplies the observed
\(\epsilon^{-1}\) factor. Gate R4 independently found that the temperature
block dominates \(\widehat c_{rr}\), while the \(G\) block partially
cancels it.

## Variable-radius linear amplitude operator

Choose \(\xi=\log(R/\ell)\), so \(X=\epsilon e^\xi\). With

\[
\widehat c_r\sim\epsilon C_r,\quad
\widehat c_{rr}\sim\epsilon^{-1}C_{rr},\quad
\widehat c_{r,\mathrm{geo}}\sim\epsilon^{-1}C_g,
\]

and

\[
A_\xi=XA_X,\quad A_{\xi\xi}=X^2A_{XX}+XA_X,\quad
e^{-2\xi}=\epsilon^2/X^2,
\]

the radial terms reduce to

\[
\boxed{
\epsilon\left[
C_{rr}A_{XX}+
\left(C_rX+\frac{C_{rr}+C_g}{X}\right)A_X
\right].
}
\]

This is the leading frozen-fold radial operator for an \(O(1)\)-\(X\)
modulation.

There is a separate amplitude-scale decision. Current finite-\(\epsilon\)
data provisionally give \(\widehat c_2=O(\epsilon^2)\), hence
\(\widehat c_2A^2=O(\epsilon^2A^2)\). After dividing by the common radial
\(O(\epsilon)\), the nonlinear term still carries \(\epsilon A^2\).
Detuning and forcing also need assigned scales. Therefore diffusion,
drift, detuning, and \(A^2\) cannot yet all be declared leading without a
consistent amplitude scaling.

The established structure is

\[
\boxed{\text{one-dimensional near-wall inner}
\longleftrightarrow
\text{two-dimensional thermo-radial outer}.}
\]

This is not yet a finite-disk solution or a complete two-dimensional
Navier--Stokes reduction.

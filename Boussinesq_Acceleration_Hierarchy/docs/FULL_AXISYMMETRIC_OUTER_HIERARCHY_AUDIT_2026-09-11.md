# Full-axisymmetric audit of the thermo-radial outer hierarchy

## Why this audit is necessary

At \(R\sim z\sim\ell/\epsilon\), the geometrical boundary-layer condition
\(z/R\ll1\) no longer holds. The \(X,Y\) outer must therefore be checked
from the full steady axisymmetric equations, not certified only by adding
radial diffusion to an already reduced boundary-layer model.

## Parent equations

Let \(\chi=\rho/\rho_\infty\), \(\Pi=p/\rho_\infty\), and subtract the
reference hydrostatic pressure when gravity is retained:

\[
\frac1R(Ru_R)_R+u_{z,z}=0,
\]

\[
\chi\left(u_Ru_{R,R}+u_zu_{R,z}-\frac{u_\phi^2}{R}\right)
=-\Pi_R+\nu\left(u_{R,RR}+\frac1Ru_{R,R}
-\frac{u_R}{R^2}+u_{R,zz}\right),
\]

\[
\chi\left(u_Ru_{\phi,R}+u_zu_{\phi,z}
+\frac{u_Ru_\phi}{R}\right)
=\nu\left(u_{\phi,RR}+\frac1Ru_{\phi,R}
-\frac{u_\phi}{R^2}+u_{\phi,zz}\right),
\]

\[
\chi(u_Ru_{z,R}+u_zu_{z,z})
=-\Pi_z+\nu\left(u_{z,RR}+\frac1Ru_{z,R}+u_{z,zz}\right)
+(1-\chi)g_z,
\]

\[
u_RT_R+u_zT_z=\alpha\left(T_{RR}+\frac1RT_R+T_{zz}\right).
\]

If gravity is not parallel to the rotation axis, axisymmetry itself is
lost.

## Distinguished substitution and pressure

Use

\[
X=\epsilon R/\ell,\quad Y=\epsilon z/\ell,
\]

\[
u_R=-\Omega\ell Ro_t\,\epsilon Xf,\quad
u_\phi=\Omega\ell\epsilon^{-1}Xq,\quad
u_z=-\Omega\ell Ro_t\,\epsilon h,
\]

where \(q=s_t+Ro_tg\), and expand

\[
\Pi=\Omega^2\ell^2\left[
\epsilon^{-2}P_0(X,Y)+\epsilon^2P_4(X,Y)+\cdots
\right].
\]

The leading radial and axial pressure balances give

\[
P_{0X}=\chi Xq^2,\qquad P_{0Y}=0.
\]

For a uniform rigid far field,

\[
P_0=\frac12\omega_t^2X^2+\text{constant},\qquad
\chi q^2=\omega_t^2.
\]

Continuity, azimuthal momentum, and energy retain the previously derived
thermo-radial equations:

\[
h_Y+2f+Xf_X=0,
\]

\[
g_{YY}+g_{XX}+\frac3Xg_X+
\chi\{Ro_t(2fg+Xfg_X+hg_Y)+Co_tf\}=0,
\]

\[
t_{YY}+t_{XX}+\frac1Xt_X+
Pr\,Ro_t(Xft_X+ht_Y)=0.
\]

These equations survive the full-equation audit, but are only a necessary
thermo-azimuthal subsystem.

## First nonzero meridional momentum order

At \(O(\Omega^2\ell\epsilon^3)\), radial and axial momentum give

\[
P_{4X}
+\chi Ro_t^2X\{f(f+Xf_X)+hf_Y\}
+Ro_t[X(f_{XX}+f_{YY})+3f_X]=0,
\]

\[
P_{4Y}
+\chi Ro_t^2(Xfh_X+hh_Y)
+Ro_t\left(h_{XX}+\frac1Xh_X+h_{YY}\right)
-\mathcal G_\epsilon(1-\chi)=0,
\]

where

\[
\mathcal G_\epsilon=\frac{g_z}{\Omega^2\ell\epsilon^3}.
\]

Although four powers below the centrifugal balance, these are the first
momentum equations governing \(f,h\). Pressure compatibility,

\[
\partial_YP_{4X}=\partial_XP_{4Y},
\]

adds a meridional-vorticity constraint. Hence the earlier four-equation
thermo-radial system is not a complete full-NS closure. This is the
principal correction from the audit.

## Gravity gate

Define the physical density contrast

\[
B_T=|\beta_\infty(T_w^{dim}-T_\infty^{dim})|
=|\gamma(T_w-1)|.
\]

Neglecting axial gravity relative to the weak meridional balance requires

\[
\boxed{\frac{gB_T}{\Omega^2\ell}\ll\epsilon^3.}
\]

This is more restrictive than \(g\ll\Omega^2R\), since meridional motion is
much weaker than the leading centrifugal balance.

## Boussinesq-validity gate

The computations use

\[
\gamma=1,\qquad T_{w,t}=1.6615874197,
\]

so

\[
\boxed{B_{T,t}=\gamma(T_{w,t}-1)=0.6615874197.}
\]

This is \(O(1)\), not an asymptotically small density change. The current
critical fold--tail must therefore be described as a formal
constant-property, linear-density organising structure until a matched
low-Mach variable-density calculation reproduces it. Changing \(\gamma\)
alone does not cure this if the equations depend only on
\(\gamma(T-1)\); it merely changes the wall-temperature coordinate.

## Finite-radius, transition and confinement parameters

Define

\[
Re_R=\frac{\Omega R^2}{\nu},\qquad
\Pi_R=\epsilon\sqrt{Re_R}=\epsilon R/\ell.
\]

Then \(\Pi_R=O(1)\) is the thermo-radial regime. If \(Re_{tr}\) is a
consistently defined heated-BEK transition threshold, the promoted region
is reached while laminar only if

\[
\boxed{\epsilon\sqrt{Re_{tr}}\gtrsim1.}
\]

No universal numerical \(Re_{tr}\) is asserted here; it depends on \(Ro\),
heating, disturbance class and Reynolds-number convention.

For a rotor--stator gap \(H_c\), define

\[
\Pi_H=\epsilon H_c/\ell.
\]

Axial confinement enters at \(\Pi_H=O(1)\). Comparing \(\Pi_H\) and
\(\Pi_R\) determines whether the opposite wall or radial geometry
regularises the thermal tail first.

## Observable validation

The minimal observables are

\[
q_w(R)=-kT_z|_{z=0},\qquad
\tau_R(R)=\mu u_{R,z}|_{z=0},\qquad
\tau_\phi(R)=\mu u_{\phi,z}|_{z=0},
\]

together with axial pumping and an edge/forcing penetration length. These
should validate the promoted scale in a weak-\(T_w(R)\) or finite-geometry
calculation; adjoint coefficients alone are not observables.

## Revised priority

1. close the full-hierarchy pressure and inner/outer boundary data;
2. build a regime map in \((\Pi_R,B_T,\Pi_H,\mathcal G)\);
3. establish whether a laminar and Boussinesq-valid window exists;
4. perform a weak radial wall-temperature forcing test;
5. only then choose between a finite-disk solve and higher-order
   matched-adjoint selection.


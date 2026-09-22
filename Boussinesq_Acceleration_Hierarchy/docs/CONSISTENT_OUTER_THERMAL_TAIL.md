# Outer thermal-tail equations of the consistent BEK model

## Scope

Derive the leading outer problem as

```text
epsilon = -Hinf -> 0+.
```

This note does not perform inner--outer matching, derive the previously fitted
`A_R,A_T`, continue any numerical branch, or add parameter slices.  Its source
model is the fully self-consistent steady residual in
`work/src/BEKConsistent.jl`.

## Governing similarity equations

Write `r=Ro`,

```text
Co = 2-r-r^2,     s = Co/2,     omega = s+r,
chi(T) = 1-gamma(T-1).
```

The steady consistent equations are

```text
H_eta + 2F = 0,

F_etaeta + chi[ r(F^2+H F_eta) - (s+rG)^2/r ]
           + omega^2/r = 0,

G_etaeta + chi[ r(2FG+H G_eta) + Co F ] = 0,

T_etaeta + Pr*Ro*H*T_eta = 0.
```

The far-field conditions are

```text
F -> 0,   G -> 1,   H -> -epsilon,   T -> 1.
```

At the fold--tail limiting point let `r=r_t+O(epsilon)`.  The `O(epsilon)`
parameter displacement is not needed in the leading outer problem.

## Temperature equation with frozen far-field velocity

First make only the far-field approximation `H=-epsilon` and set

```text
Y = epsilon eta,     T(eta)=Theta(Y)+... .
```

Then

```text
T_eta = epsilon Theta_Y,
T_etaeta = epsilon^2 Theta_YY,
```

so the temperature equation becomes

```text
Theta_YY + Pr Theta_Y = 0.
```

With `Theta(infinity)=1`,

```text
Theta(Y) = 1 + C exp(-Pr Y)
         = 1 + C exp(-Pr epsilon eta).
```

Equivalently, with `Z=Pr epsilon eta`,

```text
Theta_ZZ + Theta_Z = 0,
Theta(Z)=1+C exp(-Z).
```

Therefore the e-folding length in the original similarity coordinate is

```text
ell_T ~ 1/(Pr*(-Ro)*epsilon)  (for Ro<0 and Hinf=-epsilon).
```

The amplitude `C`, or equivalently the outer value approached as `Y->0`, is
matching data and is not determined by the outer far-field condition alone.

## Consistency check for the outer velocity

The frozen velocity is sufficient for the far-far exponential decay, but it
is not generally a uniform leading solution throughout the full `Y=O(1)`
outer layer.  If `F=0`, `G=1`, `H=-epsilon` while `Theta-1=O(1)`, the leading
radial residual is

```text
(1-chi(Theta)) omega^2/r
= gamma(Theta-1) omega^2/r.
```

It vanishes only when

```text
Theta=1,     gamma=0,     or omega=s+r=0.
```

The last special case contains `Ro=-1` (and the other algebraic endpoint
`Ro=2`), but not the present fold--tail limit near `Ro_t=-0.32`.  Hence an
`O(1)` temperature excess in the consistent outer layer induces a slow
velocity adjustment.

## Leading coupled outer problem

Use the distinguished outer ansatz

```text
Y = epsilon eta,
H = epsilon h(Y)+...,
F = epsilon^2 f(Y)+...,
G = g(Y)+...,
T = Theta(Y)+...,
r = r_t+O(epsilon).
```

Primes below denote `d/dY`; `Co,s,omega` are evaluated at `r_t`.
Continuity gives

```text
f = -h'/2.
```

The radial equation at leading order is algebraic:

```text
chi(Theta) (s+r_t g)^2 = omega^2.
```

On the branch continuous with `g->1` as `Theta->1`,

```text
s+r_t g = omega/sqrt(chi(Theta)),

g(Theta) = [omega/sqrt(chi(Theta))-s]/r_t.
```

The azimuthal equation first enters at `O(epsilon^2)`:

```text
g'' + chi(Theta)[ r_t h g' - (s+r_t g) h' ] = 0.
```

The temperature equation is

```text
Theta'' - Pr h Theta' = 0.
```

Thus the leading outer system may be written as

```text
f = -h'/2,
chi(Theta)(s+r_t g)^2 = omega^2,
g'' + chi(Theta)[r_t h g'-(s+r_t g)h'] = 0,
Theta'' - Pr h Theta' = 0,
```

with far-field conditions

```text
h -> -1,   f -> 0,   g -> 1,   Theta -> 1
as Y -> infinity.
```

Conditions as `Y->0` must be supplied by later matching and are deliberately
left unspecified here.

If `Z=Pr epsilon eta` is used instead, the last two differential equations
become

```text
Pr g_ZZ + chi(Theta)[r_t h g_Z-(s+r_t g)h_Z] = 0,
Theta_ZZ - h Theta_Z = 0.
```

## Explicit outer temperature representation

For a given outer axial-flow function `h(Y)`, integrate the temperature
equation once:

```text
Theta'(Y) = K exp(Pr integral_0^Y h(xi) dxi).
```

If `Theta(0)=Theta_m` denotes the still-unknown matching value, the solution
satisfying `Theta(infinity)=1` is

```text
Theta(Y) = 1 + (Theta_m-1)
  [ integral_Y^infinity exp(Pr integral_0^s h(xi) dxi) ds ]
  /[ integral_0^infinity exp(Pr integral_0^s h(xi) dxi) ds ].
```

When `h=-1`, this reduces exactly to

```text
Theta(Y)=1+(Theta_m-1)exp(-Pr Y).
```

More generally, because `h->-1` as `Y->infinity`, the far-far tail always has

```text
Theta-1 ~ C_infinity exp(-Pr Y)
        = C_infinity exp(-Pr epsilon eta).
```

Therefore the length `ell_T~1/(Pr epsilon)` is robust even though the full
consistent outer layer contains coupled slow variations of `G` and `H`.

## Main conclusion

- `F=O(epsilon^2)` and is zero at leading `O(1)` outer order.
- `G=1` and `H=-epsilon` are valid as `Y->infinity`, but not generally
  uniformly across the entire outer layer when `Theta-1=O(1)`.
- The far-far temperature tail is exponential with length
  `1/(Pr epsilon)`.
- The full leading outer temperature is given by the integral expression
  above and must later be closed by inner--outer matching.

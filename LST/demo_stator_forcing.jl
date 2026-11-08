#=
demo_stator_forcing.jl
======================

Illustration of the stator-side disturbance model for the rotor–stator cavity
(CavityLST.jl, rotor frame, local parallel, full gap z ∈ [0, 1]).

Stator-layer turbulence is represented by a localized fluctuating body force
(Reynolds-stress forcing) inside the stator boundary layer,
    (L0 + alpha L1 + alpha^2 L2) q = f,   f = e_k exp[-(z - z_s)^2/(2 sigma^2)],  k = r, theta, z  (z_s = 0.955, sigma = 0.02),
with no-slip on both disks. A forcing Fourier component (alpha_f, n_f, omega_f) is convected
with the local velocity of the stator layer, which in the rotor frame lies between the stator
(V = -1) and the core (V = G_c - 1): omega_f = -c beta_f, c ∈ [1 - G_c, 1]. The response is
normalised by its velocity amplitude at mid-gap, z = 0.5.
The P_N-P_{N-2} pressure formulation (`pressure=:PN2`) is required: the default P_N-P_N
collocation has a spurious zig-zag pressure mode that dominates the weak core response
(rotor-edge amplitude did not converge with N). A convergence table (N = 149 vs 249) is written.

Outputs (results/stator_forcing_demo/):
  profiles at R = 330 for several alpha_f and the three forcing directions,
  rotor-edge / core amplitude ratio versus alpha_f (penetration of the disturbance),
  stator-layer / core ratio (indicator of near-resonance with stator-layer modes).
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
using LinearAlgebra, Printf
BLAS.set_num_threads(min(16, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "stator_forcing_demo"); mkpath(OUT)

const Re_h = 1000.0
const ZS, SIG = 0.955, 0.02   # centre of the stator layer (edge at z ≈ 0.917, see printout)
const R = 330.0
const N = 199
b = cavity_baseflow(Re_h; N=300)
Gc = bary_interp(b.z, b.G, 0.5)
# rotor and stator boundary-layer edges: |G - G_c| = 0.01 (1 - G_c) on each side
zq = range(0, 1; length=4001); Gq = [bary_interp(b.z, b.G, x) for x in zq]
z_rotor_edge = zq[findfirst(g -> abs(g - Gc) < 0.01*(1 - Gc), Gq)]
z_stator_edge = zq[findlast(g -> abs(g - Gc) < 0.01*Gc, Gq)]
@printf("G_c = %.4f; rotor-layer edge z = %.3f; stator-layer edge z = %.3f; forcing centre z_s = %.2f\n", Gc, z_rotor_edge, z_stator_edge, ZS)

"""Response to a unit forcing in direction k (1 = r, 2 = theta, 3 = z) for (alpha_f, n_f, c)."""
function response(alpha_f, n_f, c, k; N=N)
    beta = n_f/R; omega = -c*beta
    op = cavity_pencil(b, R, beta, omega; N=N, pressure=:PN2)
    nn = op.n; g = exp.(-((op.z .- ZS).^2)./(2SIG^2))
    f = zeros(ComplexF64, 4nn); f[(k-1)*nn+1:k*nn] .= g
    P = op.raw[1] + alpha_f*op.raw[2] + alpha_f^2*op.raw[3]
    q = zeros(ComplexF64, 4nn); q[op.keep] .= P[op.keep, op.keep] \ f[op.keep]
    vel(zz) = sqrt(sum(abs2(bary_interp(op.z, q[(j-1)*nn+1:j*nn], zz)) for j in 1:3))
    a_core = vel(0.5)
    return (op=op, q=q ./ a_core, edge=vel(z_rotor_edge)/a_core, stator=maximum(vel(zz) for zz in range(z_stator_edge, 1, length=200))/a_core)
end

rows = NamedTuple[]
for k in 1:3, nf in (-20, -10, 0, 10), c in (0.7, 0.85, 1.0), af in vcat(0.005, 0.01:0.01:0.1, 0.15:0.05:0.6)
    r = response(af, nf, c, k)
    push!(rows, (dir=("r", "theta", "z")[k], alpha_f=af, n_f=Float64(nf), c=c, edge_over_core=r.edge, stator_over_core=r.stator))
end
open(joinpath(OUT, "penetration.csv"), "w") do f
    cols = collect(keys(rows[1])); println(f, join(string.(cols), ","))
    for r in rows; println(f, join([string(getproperty(r, c)) for c in cols], ",")); end
end
# grid convergence of the rotor-edge ratio
open(joinpath(OUT, "convergence.csv"), "w") do f
    println(f, "dir,alpha_f,n_f,c,edge_N149,edge_N249,reldiff")
    for k in 1:3, nf in (-20, 0, 10), c in (0.7, 1.0), af in (0.01, 0.1, 0.3, 0.6)
        e1 = response(af, nf, c, k; N=149).edge; e2 = response(af, nf, c, k; N=249).edge
        println(f, join((("r", "theta", "z")[k], af, nf, c, e1, e2, abs(e1 - e2)/abs(e2)), ","))
    end
end
# profiles for the figure
open(joinpath(OUT, "profiles.csv"), "w") do f
    println(f, "dir,alpha_f,z,u,v,w")
    for k in 1:3, af in (0.02, 0.1, 0.3)
        r = response(af, -10, 0.85, k); nn = r.op.n
        for (j, zz) in enumerate(r.op.z)
            println(f, join((("r", "theta", "z")[k], af, zz, abs(r.q[j]), abs(r.q[nn+j]), abs(r.q[2nn+j])), ","))
        end
    end
end
open(joinpath(OUT, "baseflow.csv"), "w") do f
    println(f, "z,F,G,H"); for j in eachindex(b.z); println(f, "$(b.z[j]),$(b.F[j]),$(b.G[j]),$(b.H[j])"); end
end
open(joinpath(OUT, "meta.txt"), "w") do f
    println(f, "Gc=$Gc"); println(f, "z_rotor_edge=$z_rotor_edge"); println(f, "z_stator_edge=$z_stator_edge"); println(f, "zs=$ZS"); println(f, "sigma=$SIG"); println(f, "R=$R")
end
println("wrote ", OUT)

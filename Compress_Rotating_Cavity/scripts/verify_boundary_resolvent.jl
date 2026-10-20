#!/usr/bin/env julia

using Printf
using Test
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))
using CompressRotatingCavity
include(joinpath(ROOT, "src", "CRC_RA.jl"))
using .CRC_RA

const N = 129
const Res = 2000.0
const Ro = -1.0
const mode = 1
const R = haskey(ENV, "RLOCAL") ? parse(Float64, ENV["RLOCAL"]) : 330.0
const beta = haskey(ENV, "BETA") ? parse(Float64, ENV["BETA"]) : 0.0503
const omega = haskey(ENV, "OMEGA") ? parse(Float64, ENV["OMEGA"]) : 8.0 / R
const alpha_seed = haskey(ENV, "ALPHA_SEED") ? parse(Float64, ENV["ALPHA_SEED"]) : 0.5
const hr = 1.0 / Res
const ls = 0.5

println("Boundary-lifted local response / residue verification")
println("  fixed base: Re=$Res, R=$R, beta=$beta, omega=$omega, N=$N")

flow = (; CompressRotatingCavity.baseflow(Res, Ro, 0.0, mode, N)..., Res)
solution = CompressRotatingCavity.normalize_modes(
    CompressRotatingCavity.eigsolve(flow, Res, Ro, mode, alpha_seed,
                                    R, omega, beta))
alpha = solution.eigval[1]
legacy_Q = CompressRotatingCavity.q_receptivity(solution,
                                                CompressRotatingCavity.cheb_quad(N))
legacy = CompressRotatingCavity.cr_receptivity(solution, flow, Res, R,
                                               hr, ls, legacy_Q)

sta = CompressRotatingCavity.CRC_STA
cof = sta.Spatial_mode_BEK1(flow.F, flow.G .- 1, flow.H, R, N,
                             flow.D, flow.D2, Res)
L0, L1, L2 = sta.assemble_mat(cof, flow.D, flow.D2, beta, omega, R)

Hx = hr * exp(-(alpha^2) / (4ls)) * sqrt(pi / ls)
uw = -(flow.D * flow.F)[1] * Hx
vw = -(flow.D * flow.G)[1] * Hx
b = wall_boundary_vector(N; uw=uw, vw=vw)
residue = spatial_residue(L0, L1, L2, N, alpha,
                         solution.eigvec[:, 1], b)

@printf("  alpha          = %.16f %+.16fi\n", real(alpha), imag(alpha))
@printf("  legacy Cr      = %.16e %+.16ei\n", real(legacy.Cr), imag(legacy.Cr))
@printf("  algebraic Cr   = %.16e %+.16ei\n", real(residue.Cr), imag(residue.Cr))
@printf("  relative diff  = %.6e\n",
        abs(residue.Cr - legacy.Cr) / max(abs(legacy.Cr), eps()))
raw_physical_abs = sqrt(Res) * abs(residue.Cr)
@printf("  raw residue after physical normal scaling |Cr| = %.16e\n",
        raw_physical_abs)
@printf("  scaled magnitude relative diff = %.6e\n",
        abs(raw_physical_abs - abs(legacy.Cr)) /
        max(abs(legacy.Cr), eps()))
if get(ENV, "CHECK_RAW", "1") == "1"
    @test isapprox(raw_physical_abs, abs(legacy.Cr); rtol=2e-3, atol=2e-8)
end

# The published Cr implementation includes the physical wall-velocity scaling
# used by the receptivity derivation.  It should be linear in the prescribed
# roughness amplitude h.  Keep this check separate from the raw collocation
# residue above: the raw P matrix uses the code's gap-scaled wall variable,
# while the paper's boundary concomitant uses the physical normal scaling.
amp0 = legacy.Cr
for scale in (0.25, 0.5, 1.0, 2.0)
    amp = CompressRotatingCavity.cr_receptivity(
        solution, flow, Res, R, hr * scale, ls, legacy_Q).Cr
    slope = amp / (hr * scale)
    @printf("  h scale=%4.2f: a_n=% .12e, a_n/h=% .12e\n",
            scale, amp, slope)
    @test isapprox(amp / (hr * scale), amp0 / hr;
                   rtol=2e-10, atol=2e-10)
end

# Verify the partitioned inhomogeneous solve at real radial wavenumbers.  The
# complex pole itself is used only for the residue; a physical Fourier forcing
# uses real kappa.
ii = interior_indices(N)
bi = boundary_indices(N)
for kappa in [0.20, 0.40, 0.60]
    P = L0 + kappa .* L1 + kappa^2 .* L2
    Hk = hr * exp(-(kappa^2) / (4ls)) * sqrt(pi / ls)
    bk = wall_boundary_vector(N;
                             uw=-(flow.D * flow.F)[1] * Hk,
                             vw=-(flow.D * flow.G)[1] * Hk)
    q = lifted_response(P, N, bk)
    interior_residual = norm(P[ii, :] * q) /
                        max(norm(P[ii, bi] * q[bi]), 1.0)
    @printf("  kappa=%.2f: BC=% .3e, interior residual=% .3e\n",
            kappa, norm(q[bi] - bk), interior_residual)
    @test norm(q[bi] - bk) < 1e-12
    @test norm(P[ii, :] * q) < 1e-9 * max(norm(P[ii, bi] * q[bi]), 1.0)
    q2 = lifted_response(P, N, 2 .* bk)
    @test norm(q2 - 2 .* q) / max(norm(q), 1.0) < 1e-12
end

# The direct/adjoint regression remains the reference for the existing model;
# this script adds the independent boundary-partition and residue check.
@test abs(solution.eigval_A[1] - alpha) < 2e-8
@test isfinite(real(residue.Cr)) && isfinite(imag(residue.Cr))
println("PASS")

#!/usr/bin/env julia

using Printf
using Test
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))
using CompressRotatingCavity
include(joinpath(ROOT, "src", "CRC_RA.jl"))
include(joinpath(ROOT, "src", "CRC_MEAN.jl"))
using .CRC_RA
using .CRC_MEAN

const N = 129
const RE = 2000.0
const R = 330.0
const AS = 0.0
const nmode = 27
const beta = nmode / R
const omega = 0.0
const alpha_seed = 0.4
const kappa = 0.4
const ls = 0.5
const h0 = 1.0 / RE

poly(L0, L1, L2, a) = L0 + a * L1 + a^2 * L2
function roughness_boundary_local(flow, N, kappa; height, ls)
    hs = height * sqrt(pi / ls) * exp(-kappa^2 / (4ls))
    return wall_boundary_vector(N;
        uw=-(flow.D * flow.F)[1] * hs,
        vw=-(flow.D * flow.G)[1] * hs)
end

println("Stationary roughness quadratic-mean forcing verification")
println("  Re_h=$RE, R=$R, n=$nmode, beta=$beta, omega=$omega, kappa=$kappa, N=$N")

flow = (; baseflow(RE, -1.0, AS, 1, N)..., Res=RE)
solution = normalize_modes(eigsolve(flow, RE, -1.0, 1, alpha_seed,
                                    R, omega, beta))
sta = CompressRotatingCavity.CRC_STA
cof = sta.Spatial_mode_BEK1(flow.F, flow.G .- 1, flow.H, R, N,
                             flow.D, flow.D2, RE)
L0, L1, L2 = sta.assemble_mat(cof, flow.D, flow.D2, beta, omega, R)
P = poly(L0, L1, L2, kappa)

function response_for_height(h)
    b = roughness_boundary_local(flow, N, kappa; height=h, ls=ls)
    return lifted_response(P, N, b)
end

function l2(v)
    return sqrt(sum(abs2, v) / length(v))
end

function hhat_for_height(h)
    return h * sqrt(pi / ls) * exp(-kappa^2 / (4ls))
end

q0 = response_for_height(h0)
f0 = mean_convective_mode(q0, flow.D, N, kappa, nmode, R)
wall20 = second_order_wall_values(flow, q0, hhat_for_height(h0),
                                  flow.D, flow.D2, N)
corr0 = mean_similarity_correction(flow, flow.D, RE,
                                   f0.radial, f0.azimuthal;
                                   wall_values=wall20)

# Frozen-R pilot: add the computed physical O(h^2) correction to the local
# similarity profiles and re-solve the existing spatial eigenproblem. This is
# deliberately diagnostic; it is not yet the full non-similar r-z base flow.
flow_corr = (; flow..., F=flow.F .+ corr0.F2,
             G=flow.G .+ corr0.G2, H=flow.H .+ corr0.H2)
solution_corr = normalize_modes(eigsolve(flow_corr, RE, -1.0, 1, real(solution.eigval[1]),
                                         R, omega, beta))
alpha0 = solution.eigval[1]
alpha_corr = solution_corr.eigval[1]
Q0 = q_receptivity(solution, cheb_quad(N))[1]
Qcorr = q_receptivity(solution_corr, cheb_quad(N))[1]
Cr0 = cr_receptivity(solution, flow, RE, R, h0, ls, [Q0]).Cr
Crcorr = cr_receptivity(solution_corr, flow_corr, RE, R, h0, ls, [Qcorr]).Cr
@printf("local frozen-R correction: |U2|=%.6e, alpha0=%.12f %+.12fi, alpha2=%.12f %+.12fi\n",
        sqrt(l2(corr0.F2)^2 + l2(corr0.G2)^2 + l2(corr0.H2)^2),
        real(alpha0), imag(alpha0), real(alpha_corr), imag(alpha_corr))
@printf("  delta alpha=%.6e %+.6ei, delta Cr=%.6e\n",
        real(alpha_corr-alpha0), imag(alpha_corr-alpha0), Crcorr-Cr0)

println("scale   ||q1||        ||a2_r||      ||a2_theta||  ||a2_z||      ||U2||")
for scale in (0.25, 0.5, 1.0, 2.0)
    q = response_for_height(scale * h0)
    f = mean_convective_mode(q, flow.D, N, kappa, nmode, R)
    wall2 = second_order_wall_values(flow, q, hhat_for_height(scale*h0),
                                     flow.D, flow.D2, N)
    corr = mean_similarity_correction(flow, flow.D, RE,
                                      f.radial, f.azimuthal;
                                      wall_values=wall2)
    u2norm = sqrt(l2(corr.F2)^2 + l2(corr.G2)^2 + l2(corr.H2)^2)
    @printf("%.2f  %.6e  %.6e  %.6e  %.6e  %.6e\n", scale, l2(q),
            l2(f.radial), l2(f.azimuthal), l2(f.axial), u2norm)
    @test isapprox(q, scale .* q0; rtol=2e-11, atol=2e-12)
    @test isapprox(f.radial, scale^2 .* f0.radial; rtol=2e-10, atol=2e-12)
    @test isapprox(f.azimuthal, scale^2 .* f0.azimuthal; rtol=2e-10, atol=2e-12)
    @test isapprox(f.axial, scale^2 .* f0.axial; rtol=2e-10, atol=2e-12)
    @test isapprox(corr.F2, scale^2 .* corr0.F2; rtol=2e-8, atol=2e-12)
    @test isapprox(corr.G2, scale^2 .* corr0.G2; rtol=2e-8, atol=2e-12)
    @test isapprox(corr.H2, scale^2 .* corr0.H2; rtol=2e-8, atol=2e-12)
    @test corr.residual < 1e-10
end

rows = ["z,mean_convective_radial,mean_convective_azimuthal,mean_convective_axial"]
for j in eachindex(flow.z)
    push!(rows, @sprintf("%.16e,%.16e,%.16e,%.16e", flow.z[j],
        f0.radial[j], f0.azimuthal[j], f0.axial[j]))
end
out = joinpath(ROOT, "benchmark_results", "mean_forcing_stationary.csv")
write(out, join(rows, "\n") * "\n")
println("Saved: $out")
println("PASS: q1, quadratic mean forcing, and frozen-similarity U2 scaling checks")

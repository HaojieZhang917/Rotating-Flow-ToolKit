#!/usr/bin/env julia

# Independent isothermal regression of the traditional thermal LST operator.
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))
using CompressRotatingCavity

const Re_h = 2000.0
const degree = 129
const R = 330.0
const beta = 0.0503
const omega = 8.0 / R
const shift = 0.5
const expected_alpha = 0.45794523797429615 - 0.01371872488304401im

println("Traditional Boussinesq LST: isothermal regression")
flow = baseflow(Re_h, -1.0, 0.0, 1, degree)
F, G, H, D, D2 = (vec(flow.F), vec(flow.G), vec(flow.H), flow.D, flow.D2)
n = degree + 1

legacy = CompressRotatingCavity.CRC_STA.Spatial_mode_BEK1(
    F, G .- 1, H, R, degree, D, D2, Re_h)
old_raw = CompressRotatingCavity.CRC_STA.assemble_mat(
    legacy, D, D2, beta, omega, R)
old_L0, old_L1, old_L2 = old_raw
old = CompressRotatingCavity.CRC_STA.boudary_condition(
    old_L0, old_L1, old_L2, degree, 1)

thermal = lst_coefficients(flow, Re_h, 0.72; radius=R)
new = spatial_matrices(thermal; beta, omega)

# At T=1, θ is decoupled.  The remaining [u,v,w,p] blocks must equal the
# paper operator exactly, before and after the wall/gauge rows are removed.
thermal_raw_L0 = thermal.L0 + beta * thermal.Lbeta + beta^2 * thermal.Lbeta2 -
                 omega * thermal.mass
thermal_raw_L1 = thermal.Lalpha
thermal_raw_L2 = thermal.Lalpha2
map_old_to_thermal = vcat(collect(1:n), collect(n+1:2n),
                          collect(2n+1:3n), collect(4n+1:5n))
old_keep = setdiff(collect(1:4n),
                   [1, n, n+1, 2n, 2n+1, 3n, 4n])
new_keep_mapped = map_old_to_thermal[old_keep]
@testset "isothermal four-field operator" begin
    @test maximum(abs, thermal_raw_L0[map_old_to_thermal, map_old_to_thermal] - old_L0) < 1e-10
    @test maximum(abs, thermal_raw_L1[map_old_to_thermal, map_old_to_thermal] - old_L1) < 1e-10
    @test maximum(abs, thermal_raw_L2[map_old_to_thermal, map_old_to_thermal] - old_L2) < 1e-10
    @test maximum(abs, thermal_raw_L0[new_keep_mapped, new_keep_mapped] - old[1]) < 1e-10
end

u = 1:n
v = n+1:2n
w = 2n+1:3n
theta = 3n+1:4n
@testset "traditional thermal columns" begin
    @test maximum(abs, thermal.L0[u, theta] + Matrix{ComplexF64}(I, n, n) / R) < 1e-14
    @test maximum(abs, thermal.L0[v, theta]) < 1e-14
    @test maximum(abs, thermal.L0[w, theta]) < 1e-14
    @test maximum(abs, thermal.L0[theta, w]) < 1e-12
end

old_val, _ = iar(PEP(collect(old)); σ=shift, neigs=1, maxit=500, tol=1e-12)
new_val, new_vec = iar(PEP(collect(new)); σ=shift, neigs=1, maxit=500, tol=1e-12)
α_old, α_new = old_val[1], new_val[1]
mode = expand_mode(thermal, new_vec[:, 1])

@printf("  legacy α = %.15f %+.15fi\n", real(α_old), imag(α_old))
@printf("  thermal α = %.15f %+.15fi\n", real(α_new), imag(α_new))
@printf("  |Δα| = %.3e; reference error = %.3e\n",
        abs(α_new - α_old), abs(α_new - expected_alpha))

@testset "isothermal eigenvalue and walls" begin
    @test abs(α_new - α_old) < 2e-8
    @test abs(α_new - expected_alpha) < 2e-8
    for field in (mode.u, mode.v, mode.w, mode.theta)
        @test abs(field[1]) == 0
        @test abs(field[end]) == 0
    end
    @test abs(mode.p[end]) == 0
end

println("PASS")

#!/usr/bin/env julia

using Printf
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))
using CompressRotatingCavity

const RE = 2000.0
const R = 330.0
const AS = 0.0
const MODE = 1
const N = parse(Int, get(ENV, "RA_N", "65"))
const n0 = parse(Int, get(ENV, "TARGET_N", "27"))
const OMEGA = parse(Float64, get(ENV, "OMEGA", "0.0"))
const ALPHA_SEED = parse(Float64, get(ENV, "ALPHA_SEED", "0.31"))
const KH = parse(Float64, get(ENV, "ROUGH_KAPPA", "0.0"))
const MLIST = [1, 2, 3, 6, 9, 12, 18, 27]

function solve_branch(flow, n::Int, seed)
    beta = n / R
    result = eigsolve(flow, RE, -1.0, MODE, seed, R, OMEGA, beta;
                      maxit=500, tol=1e-11)
    return (; n, beta, alpha=result.eigval[1], alpha_adj=result.eigval_A[1])
end

flow = (; baseflow(RE, -1.0, AS, MODE, N)..., Res=RE)
target = solve_branch(flow, n0, ALPHA_SEED)
println("Sideband detuning scan")
@printf("Re_h=%.0f, R=%.1f, N=%d, omega=%.6g, target n=%d, alpha=%.12f %+.12fi\n",
        RE, R, N, OMEGA, n0, real(target.alpha), imag(target.alpha))

rows = ["m_h,sideband,n,beta,alpha_real,alpha_imag,detuning_real,detuning_imag,detuning_abs"]
for mh in MLIST
    for sign in (-1, 1)
        ns = n0 + sign * mh
        ns == 0 && continue
        side = solve_branch(flow, ns, real(target.alpha))
        delta = side.alpha - target.alpha - KH
        push!(rows, @sprintf("%d,%+d,%d,%.12g,%.16g,%.16g,%.16g,%.16g,%.16g",
            mh, sign, ns, side.beta, real(side.alpha), imag(side.alpha),
            real(delta), imag(delta), abs(delta)))
        @printf("m_h=%2d, n=%+3d: alpha=%.9f %+.9fi, |Delta|=%.6e\n",
            mh, ns, real(side.alpha), imag(side.alpha), abs(delta))
    end
end

out = joinpath(ROOT, "benchmark_results", "sideband_detuning_stationary.csv")
write(out, join(rows, "\n") * "\n")
println("Saved: $out")

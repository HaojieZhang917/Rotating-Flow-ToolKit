#!/usr/bin/env julia

using Printf
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))
using CompressRotatingCavity

const EXPECTED_ALPHA = ComplexF64(0.45794523797429615, -0.01371872488304401)
const EXPECTED_ALPHA_A = ComplexF64(0.45794523796982584, -0.013718724883066679)
const EXPECTED_CR = 0.1392601855113145

println("Isothermal rotating-cavity benchmark")
println("  Re=2000, N=129, R=330, beta=0.0503, omega=8/330, Ts=0")

result = @timed benchmark_case()
case = result.value
alpha = case.solution.eigval[1]
alpha_A = case.solution.eigval_A[1]
Cr = case.receptivity.Cr
bc = case.flow

@printf("  elapsed   = %.3f s\n", result.time)
@printf("  base-flow BC: H(0)=%.3e, H(1)=%.3e, F(0)=%.3e, F(1)=%.3e, G(0)-1=%.3e, G(1)=%.3e\n",
        bc.H[1], bc.H[end], bc.F[1], bc.F[end], bc.G[1] - 1, bc.G[end])
@printf("  alpha     = %.16f %+.16fim\n", real(alpha), imag(alpha))
@printf("  alpha_adj = %.16f %+.16fim\n", real(alpha_A), imag(alpha_A))
@printf("  Cr        = %.16f\n", Cr)

err_alpha = abs(alpha - EXPECTED_ALPHA)
err_alpha_A = abs(alpha_A - EXPECTED_ALPHA_A)
err_Cr = abs(Cr - EXPECTED_CR)
@printf("  errors: |dalpha|=%.3e, |dalpha_adj|=%.3e, |dCr|=%.3e\n",
        err_alpha, err_alpha_A, err_Cr)

@test err_alpha < 2e-8
@test err_alpha_A < 2e-8
@test err_Cr < 2e-6
println("PASS")

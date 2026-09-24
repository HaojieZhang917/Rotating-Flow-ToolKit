#!/usr/bin/env julia

using Printf
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "CompressRotatingCavity.jl"))
using .CompressRotatingCavity

const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT_DIR = joinpath(ROOT, "benchmark_results")
const EXPECTED_TW = 1.16767648404691

function main()
    mkpath(OUTPUT_DIR)
    config = ThermalConfig(re_h=1000.0, pr=0.72, degree=100, tolerance=2e-8)
    isothermal = isothermal_profile(config, REFERENCE)
    fold, branch = trace_principal_fold(isothermal;
                                        coarse_step=0.001,
                                        fine_step=0.0002,
                                        pressure_tolerance=1e-10)
    tw_error = abs(fold.tw - EXPECTED_TW)

    branch_path = joinpath(OUTPUT_DIR, "traditional_boussinesq_Re1000_branch.csv")
    open(branch_path, "w") do io
        println(io, "pressure_gradient,Tw,residual,dTw_dPi")
        for profile in branch
            println(io, @sprintf("%.16g,%.16g,%.6e,%.6e",
                                 profile.pressure, profile.tw,
                                 profile.residual, pressure_tangent(profile)))
        end
    end

    summary_path = joinpath(OUTPUT_DIR, "traditional_boussinesq_Re1000.txt")
    open(summary_path, "w") do io
        println(io, "Traditional centrifugal Boussinesq rotor-stator benchmark")
        println(io, "Re_h=1000, Pr=0.72, degree=100, tolerance=2e-8")
        println(io, @sprintf("isothermal_pressure=%.16g", isothermal.pressure))
        println(io, @sprintf("isothermal_residual=%.6e", isothermal.residual))
        println(io, @sprintf("principal_fold_pressure=%.16g", fold.pressure))
        println(io, @sprintf("principal_fold_Tw=%.16g", fold.tw))
        println(io, @sprintf("principal_fold_residual=%.6e", fold.residual))
        println(io, @sprintf("principal_fold_dTw_dPi=%.6e", pressure_tangent(fold)))
        println(io, @sprintf("error_against_PPT=%.6e", tw_error))
        println(io, "branch_points=$(length(branch))")
    end

    @printf("Traditional centrifugal Boussinesq benchmark\n")
    @printf("  Re_h=1000, Pr=0.72, degree=100\n")
    @printf("  isothermal: Pi=%.16f, residual=%.3e\n",
            isothermal.pressure, isothermal.residual)
    @printf("  principal fold: Pi=%.16f, Tw=%.16f\n",
            fold.pressure, fold.tw)
    @printf("  residual=%.3e, dTw/dPi=%.3e, |dTw|=%.3e\n",
            fold.residual, pressure_tangent(fold), tw_error)
    @printf("  branch points=%d\n", length(branch))
    @test tw_error < 5e-8
    @test fold.residual < 1e-7
    println("PASS")
end

main()

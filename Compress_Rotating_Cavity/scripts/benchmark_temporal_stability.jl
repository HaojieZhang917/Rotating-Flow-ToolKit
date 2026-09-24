#!/usr/bin/env julia

using Printf
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
using CompressRotatingCavity

const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT_DIR = joinpath(ROOT, "benchmark_results")

function continue_to(seed, target; step=0.001)
    current = seed
    direction = sign(target - seed.pressure)
    direction == 0 && return current
    while direction * (target - current.pressure) > 1e-12
        next_pressure = current.pressure + direction * min(abs(step),
                                                           abs(target - current.pressure))
        current = fixed_pressure_profile(next_pressure, current)
    end
    return current
end

function main()
    mkpath(OUTPUT_DIR)
    config = ThermalConfig(re_h=1000.0, pr=0.72, degree=100, tolerance=2e-8)
    isothermal = isothermal_profile(config, REFERENCE)

    # Pressure values at Tw=1.160 are the three intersections of the
    # pressure-continuation curve with that temperature.
    principal_pi = 0.018230452624
    middle_pi = 0.003786919205
    upper_pi = 0.001542389588
    principal = continue_to(isothermal, principal_pi; step=0.001)
    middle = continue_to(principal, middle_pi; step=0.0002)
    upper = continue_to(middle, upper_pi; step=0.0002)
    profiles = [(1, "upper_low_Pi", upper),
                (2, "middle", middle),
                (3, "principal_isothermal", principal)]

    rows = String["branch,label,pressure_gradient,Tw,baseflow_residual,lambda_real,lambda_imag"]
    leading = Float64[]
    for (branch, label, profile) in profiles
        result = temporal_operator(profile, 120, config.re_h, config.pr)
        lambda = leading_temporal_eigenvalue(result)
        push!(leading, real(lambda))
        push!(rows, @sprintf("%d,%s,%.16g,%.16g,%.6e,%.16g,%.16g",
                            branch, label, profile.pressure, profile.tw,
                            profile.residual, real(lambda), imag(lambda)))
        for (rank, value) in enumerate(result.eigenvalues[1:min(5, end)])
            push!(rows, @sprintf("%d,%s_mode_%d,%.16g,%.16g,%.6e,%.16g,%.16g",
                                branch, label, rank, profile.pressure, profile.tw,
                                profile.residual, real(value), imag(value)))
        end
    end
    write(joinpath(OUTPUT_DIR, "traditional_boussinesq_temporal_Re1000_Tw1.160.csv"),
          join(rows, "\n") * "\n")

    summary = join([
        "Traditional Boussinesq temporal stability",
        "Re_h=1000, Pr=0.72, Tw=1.160, temporal_degree=120",
        "branch 1: upper/low-Pi",
        "branch 2: middle",
        "branch 3: principal/isothermal-connected",
        @sprintf("lambda_1=%.16g", leading[1]),
        @sprintf("lambda_2=%.16g", leading[2]),
        @sprintf("lambda_3=%.16g", leading[3]),
        "positive_real_eigenvalues=$(count(>(0.0), leading))",
    ], "\n") * "\n"
    write(joinpath(OUTPUT_DIR, "traditional_boussinesq_temporal_Re1000_Tw1.160.txt"), summary)
    print(summary)
    @test leading[1] < 0
    @test leading[2] > 0
    @test leading[3] < 0
    println("PASS")
end

main()

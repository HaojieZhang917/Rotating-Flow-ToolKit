#!/usr/bin/env julia

"""Verify how the traditional Boussinesq saddle node depends on Re_h.

For each Reynolds number, start from the isothermal state and continue in the
pressure-gradient parameter Pi toward smaller values.  The first sign change
of dTw/dPi along this branch is refined by bisection and reported as the
principal (isothermal-connected) fold.  The calculation deliberately uses Pi
as the continuation coordinate; Tw is singular at the fold.
"""

using CompressRotatingCavity
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT = joinpath(ROOT, "benchmark_results")

struct FoldResult
    re_h::Float64
    degree::Int
    status::String
    tw_isothermal::Float64
    pi_isothermal::Float64
    tw_fold::Float64
    pi_fold::Float64
    tangent_fold::Float64
    residual::Float64
    bracket_width::Float64
    steps::Int
    message::String
end

function refine_tangent_zero(lo, hi; tolerance=2e-9, max_iterations=45)
    # lo and hi bracket a zero of dTw/dPi.  The endpoint with the smaller Pi
    # is kept as lo, matching the orientation used by the continuation.
    if lo.pressure > hi.pressure
        lo, hi = hi, lo
    end
    glo = pressure_tangent(lo)
    ghi = pressure_tangent(hi)
    glo * ghi <= 0 || error("fold refinement requires a tangent sign change")
    for _ in 1:max_iterations
        hi.pressure - lo.pressure <= tolerance && break
        mid_pi = 0.5 * (lo.pressure + hi.pressure)
        seed = abs(mid_pi - lo.pressure) < abs(mid_pi - hi.pressure) ? lo : hi
        mid = fixed_pressure_profile(mid_pi, seed)
        gmid = pressure_tangent(mid)
        if glo * gmid <= 0
            hi, ghi = mid, gmid
        else
            lo, glo = mid, gmid
        end
    end
    fold_pi = 0.5 * (lo.pressure + hi.pressure)
    seed = abs(fold_pi - lo.pressure) < abs(fold_pi - hi.pressure) ? lo : hi
    fold = fixed_pressure_profile(fold_pi, seed)
    return fold, abs(hi.pressure - lo.pressure)
end

function first_principal_fold(config; pressure_step=0.0025,
                              pi_min=-0.05, tangent_tolerance=2e-9)
    iso = isothermal_profile(config, REFERENCE)
    previous = iso
    previous_tangent = pressure_tangent(previous)
    previous_tangent < 0 || error(
        "isothermal state does not start on the expected principal branch: " *
        "dTw/dPi=$(previous_tangent)",
    )
    steps = 0

    while previous.pressure > pi_min + 1e-12
        target_pi = max(pi_min, previous.pressure - pressure_step)
        current = try
            fixed_pressure_profile(target_pi, previous)
        catch exception
            # A failed large step is retried locally before declaring the
            # branch unresolved.  This keeps the production scan deterministic
            # while avoiding a hidden change of branch.
            retry_step = pressure_step / 2
            retry_step >= 1e-5 || rethrow(exception)
            target_pi = max(pi_min, previous.pressure - retry_step)
            fixed_pressure_profile(target_pi, previous)
        end
        steps += 1
        current_tangent = pressure_tangent(current)
        if previous_tangent * current_tangent <= 0
            fold, bracket_width = refine_tangent_zero(previous, current;
                                                       tolerance=tangent_tolerance)
            return iso, fold, bracket_width, steps
        end
        previous, previous_tangent = current, current_tangent
    end
    return iso, nothing, NaN, steps
end

function main()
    mkpath(OUTPUT)
    # The range includes the Re_h=1000 reference case and both lower and
    # higher Reynolds numbers.  Re_h=200 is retained to test whether the
    # principal fold leaves the positive-Pi interval.
    reynolds_numbers = (200.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0)
    degree = isempty(ARGS) ? 60 : parse(Int, ARGS[1])
    tolerance = 2e-8
    rows = FoldResult[]

    for re_h in reynolds_numbers
        println("Re_h=$(re_h)")
        config = ThermalConfig(re_h=re_h, pr=0.72, degree=degree,
                               tolerance=tolerance)
        try
            iso, fold, bracket_width, steps = first_principal_fold(config)
            if fold === nothing
                push!(rows, FoldResult(re_h, degree, "no_fold_in_scan",
                                       iso.tw, iso.pressure, NaN, NaN, NaN,
                                       iso.residual, bracket_width, steps,
                                       "No dTw/dPi sign change down to Pi=-0.05"))
                println("  no principal fold found down to Pi=-0.05")
            else
                tangent = pressure_tangent(fold)
                push!(rows, FoldResult(re_h, degree, "fold_found", iso.tw,
                                       iso.pressure, fold.tw, fold.pressure,
                                       tangent, fold.residual, bracket_width, steps,
                                       "First tangent sign change from the isothermal state"))
                println("  Tw_fold=$(fold.tw), Pi_fold=$(fold.pressure), " *
                        "residual=$(fold.residual), tangent=$(tangent)")
            end
        catch exception
            push!(rows, FoldResult(re_h, degree, "failed", NaN, NaN, NaN, NaN,
                                   NaN, NaN, NaN, 0, sprint(showerror, exception)))
            println("  FAILED: ", sprint(showerror, exception))
        end
    end

    suffix = degree == 60 ? "" : "_N$(degree)"
    path = joinpath(OUTPUT, "traditional_fold_vs_Re_h$(suffix).csv")
    open(path, "w") do io
        println(io, "Re_h,degree,status,Tw_isothermal,Pi_isothermal,Tw_fold,Pi_fold," *
                    "dTw_dPi_fold,residual,bracket_width,steps,message")
        for row in rows
            println(io, join((row.re_h, row.degree, row.status,
                              row.tw_isothermal, row.pi_isothermal,
                              row.tw_fold, row.pi_fold, row.tangent_fold,
                              row.residual, row.bracket_width, row.steps,
                              row.message), ','))
        end
    end

    summary = joinpath(OUTPUT, "traditional_fold_vs_Re_h$(suffix).txt")
    open(summary, "w") do io
        println(io, "Traditional Boussinesq principal fold versus Re_h")
        println(io, "Pr=0.72, Chebyshev degree=$(degree)")
        println(io, "The fold is the first dTw/dPi=0 point reached from Tw=1.")
        println(io, "")
        for row in rows
            println(io, join((row.re_h, row.status, row.tw_fold, row.pi_fold,
                              row.residual, row.steps, row.message), ','))
        end
    end
    println("Wrote ", path)
    println("Wrote ", summary)
end

main()

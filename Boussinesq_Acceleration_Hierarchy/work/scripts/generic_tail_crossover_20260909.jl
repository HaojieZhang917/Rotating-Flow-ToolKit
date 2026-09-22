#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using Printf

option(prefix, default) = begin
    hit = findfirst(x -> startswith(x, prefix), ARGS)
    isnothing(hit) ? default : split(ARGS[hit], "="; limit=2)[2]
end
numbers(text) = parse.(Float64, split(text, ','))
token(x) = replace(replace(@sprintf("%+.4f", x), "." => "p"), "+" => "p", "-" => "m")

const DEGREE = parse(Int, option("--degree=", "240"))
const LAMBDAS = numbers(option("--lambdas=", "-0.4,-0.3,0,0.5,1,2"))
const EPSILONS = sort(numbers(option("--eps=", "0.03,0.02,0.015,0.01")); rev=true)
const RT = parse(Float64, option("--rt=", "-0.470367416464038"))
const TT = parse(Float64, option("--tt=", "1.66158741971899"))
const MAP_A = parse(Float64, option("--a=", "2.0"))
const MAP_B = parse(Float64, option("--b=", "0.6"))
const MAP_C = parse(Float64, option("--c=", "0.5"))
const TOL = parse(Float64, option("--tol=", "1e-9"))
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTROOT = joinpath(ROOT, "work", "results", "generic_tail_crossover_20260909", "N$(DEGREE)")

function advance_h(solution, target)
    current = solution
    while current.Hinf < target - 1e-12
        trial = min(target, current.Hinf + 0.01)
        try
            current = solve_consistent_fixed_h(trial, current; tolerance=TOL)
        catch
            trial = (current.Hinf + trial)/2
            current = solve_consistent_fixed_h(trial, current; tolerance=TOL)
        end
    end
    current
end

function advance_ro(solution, target; depth=0)
    try
        solve_consistent_fixed_h_at_ro(solution.Hinf, target, solution; tolerance=TOL)
    catch
        depth >= 9 && rethrow()
        middle = (solution.Ro + target)/2
        advance_ro(advance_ro(solution, middle; depth=depth+1), target; depth=depth+1)
    end
end

function interp(eta, values, x)
    i = searchsortedlast(eta, x)
    i <= 0 && return values[1]
    i >= length(eta) && return values[end]
    t = (x-eta[i])/(eta[i+1]-eta[i])
    (1-t)*values[i] + t*values[i+1]
end

function write_profile(path, sol, epsilon, lambda)
    eta = sol.operators.eta
    H, F, G, T = eachrow(sol.fields)
    open(path, "w") do io
        println(io, "eta,H_over_eps2,F_over_eps2,G_over_eps,TminusTt_over_eps,H,F,G,T")
        for i in eachindex(eta)
            isfinite(eta[i]) || continue
            @printf(io, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                eta[i], H[i]/epsilon^2, F[i]/epsilon^2, G[i]/epsilon,
                (T[i]-TT)/epsilon, H[i], F[i], G[i], T[i])
        end
    end
end

function main()
    mkpath(OUTROOT)
    summary_path = joinpath(OUTROOT, "diagnostics.csv")
    fixed_path = joinpath(OUTROOT, "fixed_eta_scaled.csv")
    open(summary_path, "w") do sio
        println(sio, "N,Lambda,epsilon,Ro,Tw,Hinf,residual,D,D_over_eps2")
        open(fixed_path, "w") do fio
            println(fio, "N,Lambda,epsilon,Ro,eta,H_over_eps2,F_over_eps2,G_over_eps,TminusTt_over_eps")
            seed_eps0 = nothing
            for lambda in LAMBDAS
                eps0 = first(EPSILONS)
                ro0 = RT + eps0*lambda
                if isnothing(seed_eps0)
                    sol = solve_consistent_isothermal(ro0; degree=DEGREE, a=MAP_A, b=MAP_B,
                        c=MAP_C, Pr=0.72, gamma=1.0, tolerance=max(TOL,1e-8))
                    sol = advance_h(sol, -eps0)
                else
                    sol = advance_ro(seed_eps0, ro0)
                end
                seed_eps0 = sol
                for (j, epsilon) in enumerate(EPSILONS)
                    target_ro = RT + epsilon*lambda
                    if j > 1
                        sol = advance_h(sol, -epsilon)
                        sol = advance_ro(sol, target_ro)
                    end
                    H, F, G, T = eachrow(sol.fields)
                    s = sol.Co/2
                    omega = s + sol.Ro
                    D = (2-sol.Tw)*s^2 - omega^2
                    @printf(sio, "%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.5e,%.12e,%.12e\n",
                        DEGREE, lambda, epsilon, sol.Ro, sol.Tw, sol.Hinf,
                        sol.residual, D, D/epsilon^2)
                    finite = isfinite.(sol.operators.eta)
                    ef = sol.operators.eta[finite]
                    for x in (0.5, 1.0, 2.0, 5.0, 8.0)
                        @printf(fio, "%d,%.12e,%.12e,%.12e,%.4f,%.12e,%.12e,%.12e,%.12e\n",
                            DEGREE, lambda, epsilon, sol.Ro, x,
                            interp(ef,H[finite],x)/epsilon^2,
                            interp(ef,F[finite],x)/epsilon^2,
                            interp(ef,G[finite],x)/epsilon,
                            (interp(ef,T[finite],x)-TT)/epsilon)
                    end
                    pdir = joinpath(OUTROOT, "profiles", "Lambda_$(token(lambda))")
                    mkpath(pdir)
                    write_profile(joinpath(pdir, "epsilon_$(replace(@sprintf("%.3f",epsilon),"."=>"p")).csv"), sol, epsilon, lambda)
                    @printf("N=%d Lambda=%+.3f eps=%.3f Ro=%.9f residual=%.2e\n",
                        DEGREE, lambda, epsilon, sol.Ro, sol.residual)
                end
            end
        end
    end
    println("Output: $OUTROOT")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

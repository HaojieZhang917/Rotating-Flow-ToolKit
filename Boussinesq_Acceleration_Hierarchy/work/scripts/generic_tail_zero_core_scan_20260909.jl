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
token(x) = replace(@sprintf("%.4f", abs(x)), "." => "p")

const DEGREE = parse(Int, option("--degree=", "160"))
const ROS = numbers(option("--ros=", "-0.469,-0.46,-0.44,-0.40,-0.30"))
const EPSILONS = numbers(option("--eps=", "0.08,0.05,0.03,0.02,0.015,0.01"))
const MAP_A = parse(Float64, option("--a=", "2.0"))
const MAP_B = parse(Float64, option("--b=", "0.6"))
const MAP_C = parse(Float64, option("--c=", "0.5"))
const TOL = parse(Float64, option("--tol=", "1e-9"))
const TAG = option("--tag=", "")
const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TAG_SUFFIX = isempty(TAG) ? "" : "_$(TAG)"
const OUTROOT = joinpath(ROOT, "work", "results", "generic_tail_zero_core_20260909",
    "N$(DEGREE)_a$(token(MAP_A))_b$(token(MAP_B))_c$(token(MAP_C))$(TAG_SUFFIX)")

function advance_h(solution, target)
    current = solution
    step = min(0.02, max(0.0025, target-current.Hinf))
    while current.Hinf < target - 1e-12
        accepted = false
        local_step = min(step, target-current.Hinf)
        while local_step >= 2e-5
            try
                current = solve_consistent_fixed_h(current.Hinf+local_step, current;
                    tolerance=TOL)
                accepted = true
                step = min(0.02, 1.25local_step)
                break
            catch
                local_step *= 0.5
            end
        end
        accepted || error("fixed-H continuation failed before target=$target at Hinf=$(current.Hinf)")
    end
    current
end

function advance_ro(solution, target; depth=0)
    try
        return solve_consistent_fixed_h_at_ro(solution.Hinf, target, solution;
            tolerance=TOL)
    catch
        depth >= 8 && rethrow()
        midpoint = (solution.Ro + target)/2
        intermediate = advance_ro(solution, midpoint; depth=depth+1)
        return advance_ro(intermediate, target; depth=depth+1)
    end
end

function interpolate(eta, values, x)
    i = searchsortedlast(eta, x)
    i <= 0 && return values[1]
    i >= length(eta) && return values[end]
    t = (x-eta[i])/(eta[i+1]-eta[i])
    (1-t)*values[i] + t*values[i+1]
end


function inner_integrals(solution; eta_max=8.0, spacing=0.025)
    eta = solution.operators.eta
    finite = isfinite.(eta)
    ef = eta[finite]
    H, F = eachrow(solution.fields)[1:2]
    grid = collect(0.0:spacing:eta_max)
    fv = [interpolate(ef, F[finite], x) for x in grid]
    hv = [interpolate(ef, H[finite], x) for x in grid]
    trapz(values) = spacing * (sum(values) - (first(values)+last(values))/2)
    signed_h_from_f = -2trapz(fv)
    meridional_norm = sqrt(trapz(fv.^2 + hv.^2))
    return signed_h_from_f, meridional_norm
end

function write_profile(path, solution)
    eta = solution.operators.eta
    H, F, G, T = eachrow(solution.fields)
    open(path, "w") do io
        println(io, "eta,H,F,G,T,TminusTw")
        for i in eachindex(eta)
            isfinite(eta[i]) || continue
            @printf(io, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                eta[i], H[i], F[i], G[i], T[i], T[i]-solution.Tw)
        end
    end
end

function main()
    mkpath(OUTROOT)
    summary = joinpath(OUTROOT, "diagnostics.csv")
    fixed = joinpath(OUTROOT, "fixed_eta.csv")
    open(summary, "w") do sio
        println(sio, "N,a,b,c,Ro,epsilon,Tw,Hinf,residual,Fpp_wall,D_algebra,D_from_wall,identity_error,H_over_epsilon_eta8,minus2intF_over_epsilon_eta0to8,Mmer_over_epsilon_eta0to8")
        open(fixed, "w") do fio
            println(fio, "N,Ro,epsilon,eta,H,F,G,TminusTw")
            ordered_ros = sort(unique(ROS); rev=true)
            epsilon_order = sort(unique(EPSILONS); rev=true)
            base = solve_consistent_isothermal(first(ordered_ros); degree=DEGREE,
                a=MAP_A, b=MAP_B, c=MAP_C, Pr=0.72, gamma=1.0,
                tolerance=TOL)
            base = advance_h(base, -first(epsilon_order))
            seeds = Dict{Float64,ConsistentSolution}(first(ordered_ros) => base)
            previous = base
            for ro in ordered_ros[2:end]
                previous = advance_ro(previous, ro)
                seeds[ro] = previous
            end
            for ro in ordered_ros
                @printf("N=%d Ro=%.6f\n", DEGREE, ro)
                solution = seeds[ro]
                for epsilon in epsilon_order
                    solution = advance_h(solution, -epsilon)
                    H, F, G, T = eachrow(solution.fields)
                    fpp = (solution.operators.deta2 * F)[1]
                    s = solution.Co/2
                    omega = s + ro
                    chiw = 2 - solution.Tw
                    dalgebra = chiw*s^2 - omega^2
                    dwall = ro*fpp
                    h8 = interpolate(solution.operators.eta[isfinite.(solution.operators.eta)],
                        H[isfinite.(solution.operators.eta)], 8.0)/epsilon
                    signed_h, mmer = inner_integrals(solution)
                    @printf(sio, "%d,%.8f,%.8f,%.8f,%.12e,%.12e,%.12e,%.12e,%.5e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                        DEGREE, MAP_A, MAP_B, MAP_C, ro, epsilon, solution.Tw,
                        solution.Hinf, solution.residual, fpp, dalgebra, dwall,
                        dalgebra-dwall, h8, signed_h/epsilon, mmer/epsilon)
                    eta = solution.operators.eta
                    finite = isfinite.(eta)
                    ef = eta[finite]
                    for x in (1.0, 2.0, 5.0, 8.0)
                        @printf(fio, "%d,%.12e,%.12e,%.4f,%.12e,%.12e,%.12e,%.12e\n",
                            DEGREE, ro, epsilon, x,
                            interpolate(ef, H[finite], x), interpolate(ef, F[finite], x),
                            interpolate(ef, G[finite], x), interpolate(ef, T[finite], x)-solution.Tw)
                    end
                    pdir = joinpath(OUTROOT, "profiles", "Ro_m$(token(ro))")
                    mkpath(pdir)
                    write_profile(joinpath(pdir, "epsilon_$(token(epsilon)).csv"), solution)
                    @printf("  eps=%.4f Tw=%.9f D=%+.4e residual=%.2e\n",
                        epsilon, solution.Tw, dalgebra, solution.residual)
                end
            end
        end
    end
    println("Output: $OUTROOT")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

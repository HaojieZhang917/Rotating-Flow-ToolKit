#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra, DelimitedFiles, Printf

function option(prefix, default)
    for argument in ARGS
        startswith(argument, prefix) && return split(argument, "="; limit=2)[2]
    end
    default
end

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEGREE = parse(Int, option("--degree=", "240"))
const PR = 0.72
const GAMMA = 1.0
const EPSILONS = [0.02, 0.015, 0.01]
const OUTPUT_ROOT = option("--outdir=", joinpath(ROOT,"work","results",
    "corrected_manual_tangent_convergence_20260907"))
const RUN_NAMES = Dict(200=>"convergence_N200", 240=>"convergence_N240",
                       280=>"production_N280", 320=>"production_N320")

function numeric_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    names = Symbol.(split(lines[1], ','))
    [NamedTuple{Tuple(names)}(Tuple(parse.(Float64, split(line, ',')))) for line in lines[2:end]]
end

function source_fold()
    run = RUN_NAMES[DEGREE]
    directory = joinpath(ROOT, "work", "results",
        "corrected_energy_epsilon_fold_chain", run)
    rows = numeric_csv(joinpath(directory, "epsilon_fold_table.csv"))
    row = only(filter(q -> isapprox(q.epsilon, 0.02; atol=1e-14, rtol=0), rows))
    profile = readdlm(joinpath(directory, "profiles", "profile_epsilon_0p02.csv"),
                      ','; skipstart=1)
    op = BEKIsothermal._operators(DEGREE, 2.0, 0.6, 0.5)
    fields = permutedims(profile[:,2:5])
    state = BEKConsistent.pack(fields)
    residual, jacobian = BEKConsistent._residual_jacobian(
        state, row.Ro_f, row.Tw_f, GAMMA, PR, op)
    rownorm = max.([norm(view(jacobian,i,:)) for i in axes(jacobian,1)], 1e-14)
    vector = svd(jacobian ./ rownorm).V[:,end]
    vector ./= norm(vector)
    vector[length(op.x)] > 0 && (vector .*= -1)
    solution = ConsistentSolution(row.Ro_f, 2-row.Ro_f-row.Ro_f^2,
        GAMMA, PR, op, fields, row.Tw_f, row.Hinf, norm(residual,Inf), 0)
    ConsistentFold(solution, vector, norm(residual,Inf), 0)
end

function left_slope(fold)
    sol = fold.solution; op = sol.operators
    state = BEKConsistent.pack(sol.fields)
    _, jacobian = BEKConsistent._residual_jacobian(
        state, sol.Ro, sol.Tw, GAMMA, PR, op)
    rownorm = max.([norm(view(jacobian,i,:)) for i in axes(jacobian,1)], 1e-14)
    decomposition = svd(jacobian ./ rownorm)
    left = (1 ./ rownorm) .* decomposition.U[:,end]
    left ./= norm(left)
    rro = BEKConsistent._residual_ro_derivative(state, sol.Ro, sol.Tw, GAMMA, PR, op)
    rtw = BEKConsistent._tw_derivative(state, op)
    -dot(left,rro)/dot(left,rtw), decomposition.S[end], decomposition.S[end]/decomposition.S[1]
end

function analytic_slope(r)
    s = (2-r-r^2)/2; omega = s+r; sp = -0.5-r
    -2*(omega/s)*(((sp+1)*s-omega*sp)/s^2)
end

function save_mode(path, epsilon, fold)
    sol = fold.solution; n = length(sol.operators.x)
    mode = BEKConsistent.unpack(fold.nullvector,n)
    maxg = maximum(abs, view(mode,3,:))
    open(path,"w") do io
        println(io,"eta,Y,vH_scaled,vF_scaled,vG_scaled,vT_scaled")
        for i in 1:n
            @printf(io,"%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                sol.operators.eta[i], epsilon*sol.operators.eta[i],
                mode[1,i]/(maxg*epsilon), mode[2,i]/(maxg*epsilon^2),
                mode[3,i]/maxg, mode[4,i]/maxg)
        end
    end
    maximum(abs,view(mode,1,:))/maxg,
    maximum(abs,view(mode,2,:))/maxg,
    maximum(abs,view(mode,4,:))/maxg
end

function main()
    outdir = joinpath(isabspath(OUTPUT_ROOT) ? OUTPUT_ROOT : joinpath(ROOT,OUTPUT_ROOT),
        "N$(DEGREE)")
    mkpath(outdir)
    fold = source_fold()
    rows = NamedTuple[]
    for epsilon in EPSILONS
        fold = solve_consistent_fold_at_h(-epsilon, fold; tolerance=1e-9)
        tangent = fold_hinf_tangent(fold)
        adjoint, sigma, ratio = left_slope(fold)
        hratio, fratio, tratio = save_mode(joinpath(outdir,
            "mode_epsilon_$(replace(@sprintf("%.3g",epsilon),"."=>"p")).csv"),
            epsilon, fold)
        slope = tangent.Tw/tangent.Ro
        push!(rows,(;N=DEGREE,epsilon,Ro=fold.solution.Ro,Tw=fold.solution.Tw,
            dRo=tangent.Ro,dTw=tangent.Tw,slope,adjoint,
            analytic=analytic_slope(fold.solution.Ro),sigma,ratio,
            fold_residual=fold.residual,tangent_residual=tangent.residual,
            H_by_G=hratio,F_by_G=fratio,T_by_G=tratio))
        @printf("N=%d eps=%.3f Ro=%+.10f dRo=%+.9f dTw=%+.9f slope=%+.9f residual=%.2e\n",
            DEGREE,epsilon,fold.solution.Ro,tangent.Ro,tangent.Tw,slope,tangent.residual)
    end
    open(joinpath(outdir,"tangent_table.csv"),"w") do io
        println(io,"N,epsilon,Ro,Tw,dRo,dTw,dTw_dRo,adjoint_slope,analytic_slope_at_Ro,scaled_sigma_min,scaled_sigma_ratio,fold_residual,tangent_residual,Hmax_by_Gmax,Fmax_by_Gmax,Tmax_by_Gmax")
        for r in rows
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                r.N,r.epsilon,r.Ro,r.Tw,r.dRo,r.dTw,r.slope,r.adjoint,
                r.analytic,r.sigma,r.ratio,r.fold_residual,r.tangent_residual,
                r.H_by_G,r.F_by_G,r.T_by_G)
        end
    end
end

main()

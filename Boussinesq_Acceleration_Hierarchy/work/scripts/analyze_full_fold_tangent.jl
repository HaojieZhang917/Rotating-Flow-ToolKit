#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using DelimitedFiles
using Printf

function option(prefix, default)
    for argument in ARGS
        startswith(argument, prefix) && return split(argument, "="; limit=2)[2]
    end
    default
end

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_TAG = option("--results=", "corrected_energy_epsilon_fold_chain")
const DEGREE = parse(Int, option("--degree=", "240"))
const EPSILONS = parse.(Float64, split(option("--eps=", "0.02,0.01,0.005"), ','))
const CORRECT_FOLD = lowercase(option("--correct=", "false")) == "true"
const PR = 0.72
const GAMMA = 1.0
const MAP_A = 2.0
const MAP_B = 0.6
const MAP_C = 0.5

const RUN_NAMES = Dict(80=>"production_N80", 120=>"production_N120",
    160=>"production_N160", 200=>"production_N200",
    240=>"production_N240", 280=>"production_N280", 320=>"production_N320")

function numeric_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = Symbol.(split(lines[1], ','))
    [NamedTuple{Tuple(header)}(Tuple(parse.(Float64, split(line, ',')))) for line in lines[2:end]]
end

epsilon_token(epsilon) = replace(@sprintf("%.4g", epsilon), "."=>"p", "-"=>"m")

function load_fold(degree, epsilon)
    run = RUN_NAMES[degree]
    directory = joinpath(ROOT, "work", "results", RESULTS_TAG, run)
    table = numeric_csv(joinpath(directory, "epsilon_fold_table.csv"))
    row = only(filter(q -> isapprox(q.epsilon, epsilon; atol=1e-14, rtol=0), table))
    profile = readdlm(joinpath(directory, "profiles",
        "profile_epsilon_$(epsilon_token(epsilon)).csv"), ','; skipstart=1)
    op = BEKIsothermal._operators(degree, MAP_A, MAP_B, MAP_C)
    finite = isfinite.(op.eta)
    maximum(abs.(profile[finite,1] .- op.eta[finite])) < 1e-7 ||
        error("saved profile and operator grid differ")
    fields = permutedims(profile[:,2:5])
    state = BEKConsistent.pack(fields)
    residual, jacobian = BEKConsistent._residual_jacobian(
        state, row.Ro_f, row.Tw_f, GAMMA, PR, op)
    solution = ConsistentSolution(row.Ro_f, 2-row.Ro_f-row.Ro_f^2,
        GAMMA, PR, op, fields, row.Tw_f, row.Hinf, norm(residual,Inf), 0)
    rownorm = max.([norm(view(jacobian,i,:)) for i in axes(jacobian,1)], 1e-14)
    decomposition = svd(jacobian ./ rownorm)
    right = decomposition.V[:,end]
    right ./= norm(right)
    right[length(op.x)] > 0 && (right .*= -1)
    left = (1 ./ rownorm) .* decomposition.U[:,end]
    left ./= norm(left)
    fold = ConsistentFold(solution, right, norm(residual,Inf), 0)
    if CORRECT_FOLD
        try
            fold = solve_consistent_fold_at_h(-epsilon,fold;tolerance=2e-11)
        catch error
            @warn "strict fold recorrection failed; retaining saved fold" epsilon degree exception=error
        end
        solution = fold.solution
        state = BEKConsistent.pack(solution.fields)
        residual,jacobian = BEKConsistent._residual_jacobian(
            state,solution.Ro,solution.Tw,GAMMA,PR,op)
        rownorm = max.([norm(view(jacobian,i,:)) for i in axes(jacobian,1)],1e-14)
        decomposition = svd(jacobian ./ rownorm)
        left = (1 ./ rownorm) .* decomposition.U[:,end]
        left ./= norm(left)
    end
    return fold, left, decomposition.S
end

function mode_metrics(vector, n)
    parts = [view(vector, (j-1)*n+1:j*n) for j in 1:4]
    (; Hmax=maximum(abs,parts[1]), Fmax=maximum(abs,parts[2]),
       Gmax=maximum(abs,parts[3]), Tmax=maximum(abs,parts[4]),
       Hrms=norm(parts[1])/sqrt(n), Frms=norm(parts[2])/sqrt(n),
       Grms=norm(parts[3])/sqrt(n), Trms=norm(parts[4])/sqrt(n))
end

function analyze(epsilon, outdir)
    fold, left, singular = load_fold(DEGREE, epsilon)
    solution = fold.solution
    op = solution.operators
    n = length(op.x)
    m = 4n
    state = BEKConsistent.pack(solution.fields)
    residual, jacobian = BEKConsistent._residual_jacobian(
        state, solution.Ro, solution.Tw, GAMMA, PR, op)
    rro = BEKConsistent._residual_ro_derivative(
        state, solution.Ro, solution.Tw, GAMMA, PR, op)
    rtw = BEKConsistent._tw_derivative(state, op)
    solvability_slope = -dot(left,rro) / dot(left,rtw)
    tangent = fold_hinf_tangent(fold)
    tangent_slope = tangent.Tw / tangent.Ro
    # The coefficients in the historical uncorrected-energy analysis are not
    # valid after the BEK thermal advection correction.
    predicted_ro = NaN
    predicted_tw = NaN
    metrics = mode_metrics(fold.nullvector,n)

    K = BEKConsistent._jacobian_directional_derivative(
        state,solution.Ro,solution.Tw,GAMMA,PR,op,fold.nullvector)
    fold_quadratic = dot(left,K*fold.nullvector)
    tw_transversality = dot(left,rtw)

    mode = BEKConsistent.unpack(fold.nullvector,n)
    tangent_fields = BEKConsistent.unpack(tangent.field,n)
    scale = maximum(abs,view(mode,3,:))
    path = joinpath(outdir,"mode_epsilon_$(epsilon_token(epsilon)).csv")
    open(path,"w") do io
        println(io,"eta,Y,vH,vF,vG,vT,vH_by_maxG,vF_by_maxG,vG_by_maxG,vT_by_maxG,H_epsilon,F_epsilon,G_epsilon,T_epsilon")
        for i in 1:n
            @printf(io,"%.12e,%.12e",op.eta[i],epsilon*op.eta[i])
            for value in mode[:,i]; @printf(io,",%.12e",value); end
            for value in mode[:,i]./scale; @printf(io,",%.12e",value); end
            for value in tangent_fields[:,i]; @printf(io,",%.12e",value); end
            println(io)
        end
    end

    (; N=DEGREE,epsilon,Ro=solution.Ro,Tw=solution.Tw,
       base_residual=norm(residual,Inf),scaled_sigma_min=singular[end],
       scaled_sigma_ratio=singular[end]/singular[1],
       dRo=tangent.Ro,dTw=tangent.Tw,dTw_dRo=tangent_slope,
       left_slope=solvability_slope,dHinf=tangent.dHinf,
       tangent_residual=tangent.residual,
       augmented_residual=tangent.augmented_residual,
       predicted_dRo=predicted_ro,predicted_dTw=predicted_tw,
       dRo_difference=tangent.Ro-predicted_ro,
       dTw_difference=tangent.Tw-predicted_tw,
       fold_quadratic,tw_transversality,
       metrics.Hmax,metrics.Fmax,metrics.Gmax,metrics.Tmax,
       null_Hinf=fold.nullvector[n])
end

function main()
    outdir = joinpath(ROOT,"work","results","corrected_energy_full_fold_tangent","N$(DEGREE)")
    mkpath(outdir)
    rows = NamedTuple[]
    for epsilon in EPSILONS
        @printf("N=%d epsilon=%.6g: loading fold and solving augmented tangent ...\n",DEGREE,epsilon)
        row = analyze(epsilon,outdir)
        push!(rows,row)
        @printf("  dRo/deps=%+.10f dTw/deps=%+.10f ratio=%+.10f left=%+.10f residual=%.2e\n",
                row.dRo,row.dTw,row.dTw_dRo,row.left_slope,row.tangent_residual)
    end
    open(joinpath(outdir,"tangent_table.csv"),"w") do io
        names = propertynames(rows[1])
        println(io,join(string.(names),','))
        for row in rows
            for (j,name) in enumerate(names)
                j > 1 && print(io,',')
                value = getproperty(row,name)
                if value isa Integer
                    print(io,value)
                else
                    @printf(io,"%.12e",value)
                end
            end
            println(io)
        end
    end
    open(joinpath(outdir,"scope.txt"),"w") do io
        println(io,"Inputs are saved finite-epsilon folds; no new continuation or epsilon point.")
        println(io,"The tangent is obtained by differentiating the complete augmented fold system.")
        println(io,"epsilon=-Hinf fixes dHinf/depsilon=-1.")
    end
    println("Outputs: ",outdir)
end

main()

#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using Printf
using DelimitedFiles
using Plots

function option(prefix,default)
    for argument in ARGS
        startswith(argument,prefix) &&
            return split(argument,"=";limit=2)[2]
    end
    default
end

const DEGREE = parse(Int,option("--degree=","120"))
const TAG = option("--tag=","production")
const TOLERANCE = parse(Float64,option("--tol=","1e-9"))
const MAP_A = parse(Float64,option("--map-a=","2.0"))
const MAP_B = parse(Float64,option("--map-b=","0.6"))
const MAP_C = parse(Float64,option("--map-c=","0.5"))
const PRANDTL = 0.72
const GAMMA = 1.0
const EPSILONS = parse.(Float64,split(option("--eps=","0.1,0.05,0.02,0.01,0.005"),','))

function approach_tw(ro,targets)
    solution = solve_consistent_isothermal(ro;degree=DEGREE,a=MAP_A,
        b=MAP_B,c=MAP_C,Pr=PRANDTL,gamma=GAMMA,
        tolerance=TOLERANCE)
    for tw in targets
        solution = solve_consistent_fixed_tw(tw,solution;
                                               tolerance=TOLERANCE)
    end
    solution
end

function initial_fold()
    solution = approach_tw(-0.5,Float64[])
    target = -0.1
    nominal_step = 0.01
    minimum_step = 1e-5
    while solution.Hinf < target-1e-12
        proposed = min(nominal_step,target-solution.Hinf)
        local_step = proposed
        accepted = false
        while local_step >= minimum_step
            hnext = min(target,solution.Hinf+local_step)
            try
                solution = solve_consistent_fixed_h(hnext,solution;
                                                     tolerance=TOLERANCE)
                accepted = true
                nominal_step = min(0.02,1.25local_step)
                break
            catch
                local_step *= 0.5
            end
        end
        accepted || error(@sprintf(
            "fixed-H seed path failed before Hinf=%.8g; last Hinf=%.8g",
            target,solution.Hinf))
    end
    # At high degree the augmented fold Newton solve needs a seed much closer
    # to the fixed-Ro turning point than Hinf=-0.1.  Bracket it cheaply along
    # the regular fixed-H branch and select the largest Tw state.
    candidates = [solution]
    for h in -0.098:0.002:-0.07
        solution = solve_consistent_fixed_h(h,solution;
                                             tolerance=TOLERANCE)
        push!(candidates,solution)
    end
    seed = candidates[argmax(getfield.(candidates,:Tw))]
    solve_consistent_fold(seed;tolerance=TOLERANCE)
end

function advance_to_h(fold,target)
    current = fold
    minimum_step = 2e-6
    while abs(current.solution.Hinf-target) > 1e-12
        gap = target-current.solution.Hinf
        direction = sign(gap)
        proposed = min(abs(gap),0.02,max(5e-4,0.2abs(current.solution.Hinf)))
        accepted = false
        local_step = proposed
        while local_step >= minimum_step
            hnext = current.solution.Hinf + direction*local_step
            direction > 0 && (hnext = min(target,hnext))
            direction < 0 && (hnext = max(target,hnext))
            try
                candidate = solve_consistent_fold_at_h(hnext,current;
                    tolerance=TOLERANCE)
                current = candidate
                accepted = true
                break
            catch
                local_step *= 0.5
            end
        end
        accepted || error(@sprintf(
            "fixed-H fold corrector failed before Hinf=%.8g; last Hinf=%.8g",
            target,current.solution.Hinf))
    end
    current
end

function diagnostics(fold)
    solution = fold.solution
    state = BEKConsistent.pack(solution.fields)
    _,jacobian = BEKConsistent._residual_jacobian(state,solution.Ro,
        solution.Tw,solution.gamma,solution.Pr,solution.operators)
    rownorm = max.([norm(view(jacobian,row,:))
                    for row in axes(jacobian,1)],1e-14)
    decomposition = svd(jacobian ./ rownorm)
    left = (1 ./ rownorm).*decomposition.U[:,end]
    left ./= norm(left)
    vector = fold.nullvector
    K = BEKConsistent._jacobian_directional_derivative(state,solution.Ro,
        solution.Tw,solution.gamma,solution.Pr,solution.operators,vector)
    rtw = BEKConsistent._tw_derivative(state,solution.operators)
    (sigma_min=decomposition.S[end],
     sigma_ratio=decomposition.S[end]/decomposition.S[1],
     fold_quadratic=dot(left,K*vector),
     tw_transversality=dot(left,rtw),
     null_hinf=vector[length(solution.operators.x)])
end

function epsilon_token(epsilon)
    replace(@sprintf("%.4g",epsilon),"."=>"p","-"=>"m")
end

function write_profile(path,fold)
    solution = fold.solution
    data = hcat(solution.operators.eta,solution.fields')
    open(path,"w") do io
        println(io,"eta,H,F,G,T")
        writedlm(io,data,',')
    end
end

function main()
    output = joinpath(@__DIR__,"..","results",
        "corrected_energy_epsilon_fold_chain","$(TAG)_N$(DEGREE)")
    profiles = joinpath(output,"profiles")
    mkpath(profiles)
    println("Correcting the Ro=-0.5 fold seed ...")
    fold = initial_fold()
    rows = NamedTuple[]
    for epsilon in EPSILONS
        target = -epsilon
        @printf("epsilon=%.6g (Hinf=%+.6g) ...\n",epsilon,target)
        fold = advance_to_h(fold,target)
        d = diagnostics(fold)
        solution = fold.solution
        row = (epsilon=epsilon,Ro=solution.Ro,Tw=solution.Tw,
            Hinf=solution.Hinf,
            thermal_length=1/(-solution.Pr*solution.Ro*epsilon),
            sigma_min=d.sigma_min,sigma_ratio=d.sigma_ratio,
            fold_quadratic=d.fold_quadratic,
            tw_transversality=d.tw_transversality,
            null_hinf=d.null_hinf,residual=fold.residual,
            iterations=fold.iterations)
        push!(rows,row)
        write_profile(joinpath(profiles,
            "profile_epsilon_$(epsilon_token(epsilon)).csv"),fold)
        @printf("  Ro=%+.12f Tw=%.12f sigma=%.3e residual=%.3e\n",
                row.Ro,row.Tw,row.sigma_min,row.residual)
    end

    open(joinpath(output,"epsilon_fold_table.csv"),"w") do io
        println(io,"epsilon,Ro_f,Tw_f,Hinf,thermal_length,scaled_sigma_min,scaled_sigma_ratio,fold_quadratic,tw_transversality,nullvector_Hinf,residual,iterations")
        for row in rows
            @printf(io,"%.12e,%.12f,%.12f,%.12f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.5e,%d\n",
                    row.epsilon,row.Ro,row.Tw,row.Hinf,row.thermal_length,
                    row.sigma_min,row.sigma_ratio,row.fold_quadratic,
                    row.tw_transversality,row.null_hinf,row.residual,
                    row.iterations)
        end
    end

    open(joinpath(output,"run_summary.txt"),"w") do io
        @printf(io,"degree=%d\na=%.6f\nb=%.6f\nc=%.6f\nPr=%.6f\ngamma=%.6f\n",
                DEGREE,MAP_A,MAP_B,MAP_C,PRANDTL,GAMMA)
        @printf(io,"points=%d\n",length(rows))
        if length(rows) >= 2
            @printf(io,"last_delta_Ro=%.12e\n",rows[end].Ro-rows[end-1].Ro)
            @printf(io,"last_delta_Tw=%.12e\n",rows[end].Tw-rows[end-1].Tw)
        end
    end

    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8,legend=false)
    p1 = plot(getfield.(rows,:epsilon),getfield.(rows,:Ro);
              xscale=:log10,xflip=true,xlabel="epsilon=-Hinf",ylabel="Ro_f")
    p2 = plot(getfield.(rows,:epsilon),getfield.(rows,:Tw);
              xscale=:log10,xflip=true,xlabel="epsilon=-Hinf",ylabel="T_w,f")
    savefig(plot(p1,p2;layout=(1,2),size=(900,360)),
            joinpath(output,"epsilon_fold_chain.png"))
    println("Outputs: $output")
end

main()

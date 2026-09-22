#!/usr/bin/env julia

using DelimitedFiles
using LinearAlgebra
using Printf
using Statistics
using Plots

const ROOT = joinpath(@__DIR__,"..","results",
                      "consistent_epsilon_fold_chain")
const RUNS = [(80,"convergence_N80"),(120,"production_N120"),
              (160,"convergence_N160"),(200,"convergence_N200"),
              (240,"convergence_N240"),(280,"convergence_N280"),
              (320,"convergence_N320")]
const HIGH_DEGREES = Set([240,280,320])

function numeric_csv(path)
    lines = readlines(path)
    header = Symbol.(split(first(lines),','))
    rows = NamedTuple[]
    for line in lines[2:end]
        values = parse.(Float64,split(line,','))
        push!(rows,NamedTuple{Tuple(header)}(Tuple(values)))
    end
    rows
end

function profile_path(directory,epsilon)
    token = replace(@sprintf("%.4g",epsilon),"."=>"p","-"=>"m")
    joinpath(ROOT,directory,"profiles","profile_epsilon_$(token).csv")
end

function main()
    allrows = NamedTuple[]
    bydegree = Dict{Int,Vector{NamedTuple}}()
    for (degree,directory) in RUNS
        rows = numeric_csv(joinpath(ROOT,directory,"epsilon_fold_table.csv"))
        bydegree[degree] = rows
        for row in rows
            push!(allrows,(N=degree,epsilon=row.epsilon,Ro=row.Ro_f,
                           Tw=row.Tw_f,residual=row.residual,
                           thermal_length=row.thermal_length))
        end
    end
    open(joinpath(ROOT,"epsilon_fold_convergence.csv"),"w") do io
        println(io,"N,epsilon,Ro_f,Tw_f,residual,thermal_length")
        for row in allrows
            @printf(io,"%d,%.12e,%.12f,%.12f,%.5e,%.12e\n",
                    row.N,row.epsilon,row.Ro,row.Tw,row.residual,
                    row.thermal_length)
        end
    end

    epsilons = getfield.(bydegree[320],:epsilon)
    recommended = NamedTuple[]
    for epsilon in epsilons
        selected = [row for row in allrows
                    if row.N in HIGH_DEGREES && row.epsilon == epsilon]
        ros = getfield.(selected,:Ro)
        tws = getfield.(selected,:Tw)
        push!(recommended,(epsilon=epsilon,Ro=mean(ros),Tw=mean(tws),
            Ro_uncertainty=(maximum(ros)-minimum(ros))/2,
            Tw_uncertainty=(maximum(tws)-minimum(tws))/2,
            N_min=minimum(getfield.(selected,:N)),
            N_max=maximum(getfield.(selected,:N))))
    end
    open(joinpath(ROOT,"recommended_epsilon_fold_table.csv"),"w") do io
        println(io,"epsilon,Ro_f,Tw_f,Ro_half_range,Tw_half_range,N_consensus")
        for row in recommended
            @printf(io,"%.12e,%.12f,%.12f,%.12e,%.12e,%d-%d\n",
                    row.epsilon,row.Ro,row.Tw,row.Ro_uncertainty,
                    row.Tw_uncertainty,row.N_min,row.N_max)
        end
    end

    # This is only a finite-limit diagnostic: fit the last three computed
    # points linearly in epsilon and report the intercept without assigning
    # an asymptotic order or scaling law.
    tail = recommended[end-2:end]
    design = hcat(ones(length(tail)),getfield.(tail,:epsilon))
    rofit = design\getfield.(tail,:Ro)
    twfit = design\getfield.(tail,:Tw)
    open(joinpath(ROOT,"finite_limit_diagnostic.txt"),"w") do io
        println(io,"method=linear_intercept_of_last_three_computed_epsilon_points")
        println(io,"scope=finite-limit diagnostic only; not an asymptotic derivation")
        @printf(io,"Ro_intercept=%.12f\n",rofit[1])
        @printf(io,"Tw_intercept=%.12f\n",twfit[1])
        @printf(io,"last_step_delta_Ro=%.12e\n",tail[end].Ro-tail[end-1].Ro)
        @printf(io,"last_step_delta_Tw=%.12e\n",tail[end].Tw-tail[end-1].Tw)
    end

    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8,legendfontsize=8)
    eps = getfield.(recommended,:epsilon)
    p1 = plot(eps,getfield.(recommended,:Ro);
        yerror=getfield.(recommended,:Ro_uncertainty),xscale=:log10,
        xflip=true,marker=:circle,label="N=240--320 consensus",
        xlabel="epsilon=-Hinf",ylabel="Ro_f")
    p2 = plot(eps,getfield.(recommended,:Tw);
        yerror=getfield.(recommended,:Tw_uncertainty),xscale=:log10,
        xflip=true,marker=:circle,label="N=240--320 consensus",
        xlabel="epsilon=-Hinf",ylabel="T_w,f")
    savefig(plot(p1,p2;layout=(1,2),size=(920,360)),
            joinpath(ROOT,"recommended_epsilon_fold_chain.png"))

    pH = plot(;xlabel="eta",ylabel="H",xlim=(0,30))
    pT = plot(;xlabel="1+eta (log scale)",
              ylabel="(T-1)/(Tw-1)",xscale=:log10,xlim=(1,2000))
    colors = cgrad(:viridis,length(epsilons),categorical=true)
    for (index,epsilon) in enumerate(epsilons)
        data = readdlm(profile_path("convergence_N320",epsilon),',';
                       skipstart=1)
        eta,H,T = data[:,1],data[:,2],data[:,5]
        finite = isfinite.(eta)
        label = @sprintf("epsilon=%.3g",epsilon)
        plot!(pH,eta[finite],H[finite];color=colors[index],label=label)
        Tw = first(T)
        theta = (T.-1)./(Tw-1)
        plot!(pT,1 .+ eta[finite],theta[finite];
              color=colors[index],label=label)
    end
    savefig(plot(pH,pT;layout=(1,2),size=(960,380)),
            joinpath(ROOT,"epsilon_fold_profiles.png"))

    @printf("finite-limit diagnostic: Ro=%.9f Tw=%.9f\n",
            rofit[1],twfit[1])
    for row in recommended
        @printf("epsilon=%7.4f Ro=%+.9f +/- %.1e Tw=%.9f +/- %.1e\n",
                row.epsilon,row.Ro,row.Ro_uncertainty,row.Tw,
                row.Tw_uncertainty)
    end
end

main()

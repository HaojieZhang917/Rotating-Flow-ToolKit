#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Plots

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const INROOT = joinpath(ROOT, "work", "results", "generic_tail_crossover_20260909")
const OUT = joinpath(INROOT, "analysis")
const RT = -0.470367416464038

function rows(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(first(lines), ',')
    [Dict(header .=> split(line, ',')) for line in lines[2:end]]
end
num(r, key) = parse(Float64, r[key])

function grouped_fixed(degree)
    data = rows(joinpath(INROOT, "N$(degree)", "fixed_eta_scaled.csv"))
    groups = Dict{Tuple{Float64,Float64},Vector{Dict{SubString{String},SubString{String}}}}()
    for r in data
        push!(get!(groups, (num(r,"Lambda"),num(r,"epsilon")), []), r)
    end
    groups
end

function collapse_table(degree)
    groups = grouped_fixed(degree)
    lambdas = sort(unique(first.(keys(groups))))
    fields = ("H_over_eps2", "F_over_eps2", "G_over_eps", "TminusTt_over_eps")
    output = NamedTuple[]
    for lambda in lambdas
        reference = groups[(lambda,0.01)]
        for epsilon in (0.03,0.02,0.015)
            current = groups[(lambda,epsilon)]
            for field in fields
                y = [num(r,field) for r in current]
                yr = [num(r,field) for r in reference]
                rel = norm(y-yr)/max(norm(yr),1e-14)
                push!(output,(;degree,lambda,epsilon,field,relative_l2=rel))
            end
        end
    end
    output
end

function grid_table()
    a, b = grouped_fixed(240), grouped_fixed(320)
    fields = ("H_over_eps2", "F_over_eps2", "G_over_eps", "TminusTt_over_eps")
    output = NamedTuple[]
    for key in sort(collect(keys(a)))
        for field in fields
            y = [num(r,field) for r in a[key]]
            yr = [num(r,field) for r in b[key]]
            rel = norm(y-yr)/max(norm(yr),1e-14)
            push!(output,(;lambda=key[1],epsilon=key[2],field,relative_l2=rel))
        end
    end
    output
end

function fold_table()
    paths = Dict(
        280 => joinpath(ROOT,"work","results","corrected_energy_epsilon_fold_chain","production_N280","epsilon_fold_table.csv"),
        320 => joinpath(ROOT,"work","results","corrected_energy_epsilon_fold_chain","production_N320","epsilon_fold_table.csv"))
    output = NamedTuple[]
    for (degree,path) in paths, r in rows(path)
        epsilon=num(r,"epsilon"); ro=num(r,"Ro_f")
        push!(output,(;degree,epsilon,ro,lambda=(ro-RT)/epsilon))
    end
    sort(output,by=x->(x.epsilon,x.degree);rev=true)
end

function defect_fits()
    output=NamedTuple[]
    for degree in (240,320)
        data=rows(joinpath(INROOT,"N$(degree)","diagnostics.csv"))
        for epsilon in sort(unique(num(r,"epsilon") for r in data);rev=true)
            d=filter(r->abs(num(r,"epsilon")-epsilon)<1e-12,data)
            x=[num(r,"Lambda") for r in d]
            y=[num(r,"D_over_eps2") for r in d]
            X=hcat(ones(length(x)),x)
            c=X\y
            rms=norm(y-X*c)/sqrt(length(y))
            push!(output,(;degree,epsilon,intercept=c[1],slope=c[2],rms))
        end
    end
    output
end

function main()
    mkpath(OUT)
    collapses=vcat(collapse_table(240),collapse_table(320))
    open(joinpath(OUT,"fixed_lambda_collapse_errors.csv"),"w") do io
        println(io,"N,Lambda,epsilon,field,relative_L2_vs_epsilon_0p01")
        for x in collapses
            @printf(io,"%d,%.6f,%.6f,%s,%.12e\n",x.degree,x.lambda,x.epsilon,x.field,x.relative_l2)
        end
    end
    grids=grid_table()
    open(joinpath(OUT,"grid_convergence.csv"),"w") do io
        println(io,"Lambda,epsilon,field,relative_L2_N240_vs_N320")
        for x in grids
            @printf(io,"%.6f,%.6f,%s,%.12e\n",x.lambda,x.epsilon,x.field,x.relative_l2)
        end
    end
    folds=fold_table()
    open(joinpath(OUT,"fold_crossover_coordinate.csv"),"w") do io
        println(io,"N,epsilon,Ro_f,Lambda_f")
        for x in folds
            @printf(io,"%d,%.12e,%.12e,%.12e\n",x.degree,x.epsilon,x.ro,x.lambda)
        end
    end
    fits=defect_fits()
    open(joinpath(OUT,"defect_affine_crossover_fit.csv"),"w") do io
        println(io,"N,epsilon,intercept,slope_in_Lambda,rms_residual")
        for x in fits
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e\n",
                x.degree,x.epsilon,x.intercept,x.slope,x.rms)
        end
    end

    p=plot(;xlabel="epsilon",ylabel="Lambda_f = (Ro_f-Ro_t)/epsilon",
        xscale=:log10,framestyle=:box,gridalpha=.25,size=(760,470))
    for degree in (280,320)
        d=sort(filter(x->x.degree==degree,folds),by=x->x.epsilon)
        plot!(getfield.(d,:epsilon),getfield.(d,:lambda),marker=:circle,label="N=$degree")
    end
    savefig(joinpath(OUT,"fold_enters_crossover_layer.png"))

    fixed=rows(joinpath(INROOT,"N320","fixed_eta_scaled.csv"))
    p2=plot(;xlabel="Lambda",ylabel="scaled inner field at eta=5",
        framestyle=:box,gridalpha=.25,size=(800,500))
    for (field,label) in (("H_over_eps2","H/eps^2"),("F_over_eps2","F/eps^2"))
        for epsilon in (0.03,0.02,0.015,0.01)
            d=sort(filter(r->abs(num(r,"eta")-5)<1e-9 && abs(num(r,"epsilon")-epsilon)<1e-12,fixed),by=r->num(r,"Lambda"))
            plot!([num(r,"Lambda") for r in d],[num(r,field) for r in d],marker=:circle,
                label="$label, eps=$epsilon")
        end
    end
    savefig(joinpath(OUT,"crossover_scaled_fields_eta5.png"))
    println("Output: $OUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

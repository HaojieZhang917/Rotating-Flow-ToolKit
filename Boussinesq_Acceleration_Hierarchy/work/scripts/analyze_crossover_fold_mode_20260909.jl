#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Plots

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const SRC=joinpath(ROOT,"work","results","presentation_evidence_audit_20260908","tangent")
const OUT=joinpath(ROOT,"work","results","crossover_fold_solvability_20260909")

function rows(path)
    lines=filter(!isempty,strip.(readlines(path)))
    header=split(first(lines),',')
    [Dict(header .=> split(line,',')) for line in lines[2:end]]
end
num(r,k)=parse(Float64,r[k])

function interp(x,y,target)
    keep=isfinite.(x) .& isfinite.(y)
    xx=x[keep]; yy=y[keep]
    i=searchsortedlast(xx,target)
    i<=0 && return yy[1]
    i>=length(xx) && return yy[end]
    t=(target-xx[i])/(xx[i+1]-xx[i])
    (1-t)*yy[i]+t*yy[i+1]
end

mode_token(eps)=eps==0.02 ? "0p02" : eps==0.015 ? "0p015" : "0p01"

function main()
    mkpath(OUT)
    data=NamedTuple[]
    for N in (200,240,280,320)
        path=joinpath(SRC,"N$(N)","tangent_table.csv")
        isfile(path) || continue
        for r in rows(path)
            eps=num(r,"epsilon")
            push!(data,(;N,epsilon=eps,Ro=num(r,"Ro"),Tw=num(r,"Tw"),
                dRo=num(r,"dRo"),dTw=num(r,"dTw"),
                Hscaled=num(r,"Hmax_by_Gmax")/eps,
                Fscaled=num(r,"Fmax_by_Gmax")/eps^2,
                Tscaled=num(r,"Tmax_by_Gmax"),
                residual=num(r,"tangent_residual")))
        end
    end
    sort!(data,by=x->(x.N,-x.epsilon))
    open(joinpath(OUT,"null_mode_scaling.csv"),"w") do io
        println(io,"N,epsilon,Ro,Tw,dRo_depsilon,dTw_depsilon,Hmax_over_Gmax_epsilon,Fmax_over_Gmax_epsilon2,Tmax_over_Gmax,tangent_residual")
        for x in data
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                x.N,x.epsilon,x.Ro,x.Tw,x.dRo,x.dTw,x.Hscaled,x.Fscaled,x.Tscaled,x.residual)
        end
    end
    open(joinpath(OUT,"N320_linear_extrapolation.csv"),"w") do io
        println(io,"quantity,limit_at_epsilon_zero,linear_coefficient,rms_residual")
        d=sort(filter(x->x.N==320,data),by=x->x.epsilon)
        X=hcat(ones(length(d)),getfield.(d,:epsilon))
        for q in (:dRo,:dTw,:Hscaled,:Fscaled,:Tscaled)
            y=getfield.(d,q); c=X\y; rms=norm(y-X*c)/sqrt(length(y))
            @printf(io,"%s,%.12e,%.12e,%.12e\n",String(q),c[1],c[2],rms)
        end
    end
    open(joinpath(OUT,"inner_fixed_eta_scaling.csv"),"w") do io
        println(io,"N,epsilon,eta,vH_over_Geps2,vF_over_Geps2,vG_over_Geps,vT_over_Geps")
        for N in (240,280,320), eps in (0.02,0.015,0.01)
            m=rows(joinpath(SRC,"N$(N)","mode_epsilon_$(mode_token(eps)).csv"))
            eta=[num(r,"eta") for r in m]
            hs=[num(r,"vH_scaled") for r in m]
            fs=[num(r,"vF_scaled") for r in m]
            gs=[num(r,"vG_scaled") for r in m]
            ts=[num(r,"vT_scaled") for r in m]
            for x in (1.0,2.0,5.0)
                @printf(io,"%d,%.12e,%.4f,%.12e,%.12e,%.12e,%.12e\n",
                    N,eps,x,interp(eta,hs,x)/eps,interp(eta,fs,x),
                    interp(eta,gs,x)/eps,interp(eta,ts,x)/eps)
            end
        end
    end
    open(joinpath(OUT,"outer_fixed_Y_scaling.csv"),"w") do io
        println(io,"N,epsilon,Y,vH_over_Geps,vF_over_Geps2,vG_over_G,vT_over_G")
        for N in (240,280,320), eps in (0.02,0.015,0.01)
            m=rows(joinpath(SRC,"N$(N)","mode_epsilon_$(mode_token(eps)).csv"))
            Y=[num(r,"Y") for r in m]
            hs=[num(r,"vH_scaled") for r in m]
            fs=[num(r,"vF_scaled") for r in m]
            gs=[num(r,"vG_scaled") for r in m]
            ts=[num(r,"vT_scaled") for r in m]
            for x in (0.1,0.25,0.5,1.0)
                @printf(io,"%d,%.12e,%.4f,%.12e,%.12e,%.12e,%.12e\n",
                    N,eps,x,interp(Y,hs,x),interp(Y,fs,x),
                    interp(Y,gs,x),interp(Y,ts,x))
            end
        end
    end
    p=plot(;xlabel="epsilon",ylabel="scaled fold-null amplitude",
        framestyle=:box,gridalpha=.25,size=(800,500))
    for N in (240,280,320)
        d=sort(filter(x->x.N==N,data),by=x->x.epsilon)
        plot!(getfield.(d,:epsilon),getfield.(d,:Hscaled),marker=:circle,label="H/(G eps), N=$N")
        plot!(getfield.(d,:epsilon),getfield.(d,:Fscaled),marker=:square,label="F/(G eps^2), N=$N")
    end
    savefig(joinpath(OUT,"null_mode_scaled_amplitudes.png"))
    println("Output: $OUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Statistics
using Plots

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const TANGENT_ROOT=joinpath(ROOT,"work","results","full_fold_tangent_analysis")
const STUDY_ROOT=joinpath(ROOT,"work","results","tangent_convergence_mode_collapse")
const DEGREES=[200,240,280,320,360]
const EPSILONS=[0.02,0.015,0.01]

function dict_csv(path)
    lines=filter(!isempty,strip.(readlines(path)))
    header=split(lines[1],',')
    [Dict(header .=> parse.(Float64,split(line,','))) for line in lines[2:end]]
end

function findrow(rows,epsilon)
    only(filter(r->isapprox(r["epsilon"],epsilon;atol=1e-14,rtol=0),rows))
end

function tangent_rows()
    output=NamedTuple[]
    for N in DEGREES
        if N<360
            existing=dict_csv(joinpath(TANGENT_ROOT,"N$(N)","tangent_table.csv"))
            middle=dict_csv(joinpath(STUDY_ROOT,"N$(N)","fold_and_tangent.csv"))
            for epsilon in EPSILONS
                if epsilon==0.015
                    r=findrow(middle,epsilon)
                    push!(output,(;N,epsilon,dRo=r["dRo_depsilon"],dTw=r["dTw_depsilon"],
                        ratio=r["dTw_dRo"],residual=r["tangent_residual"]))
                else
                    r=findrow(existing,epsilon)
                    push!(output,(;N,epsilon,dRo=r["dRo"],dTw=r["dTw"],
                        ratio=r["dTw_dRo"],residual=r["tangent_residual"]))
                end
            end
        else
            rows=dict_csv(joinpath(STUDY_ROOT,"N360","tangent_table.csv"))
            for epsilon in EPSILONS
                r=findrow(rows,epsilon)
                push!(output,(;N,epsilon,dRo=r["dRo_depsilon"],dTw=r["dTw_depsilon"],
                    ratio=r["dTw_dRo"],residual=r["tangent_residual"]))
            end
        end
    end
    output
end

function linear_fit(x,y)
    hcat(ones(length(x)),x)\y
end

function interpolation(x,y,xnew)
    output=similar(xnew)
    for i in eachindex(xnew)
        q=xnew[i]
        j=clamp(searchsortedlast(x,q),1,length(x)-1)
        weight=(q-x[j])/(x[j+1]-x[j])
        output[i]=(1-weight)*y[j]+weight*y[j+1]
    end
    output
end

function mode_collapse()
    mode_rows=Dict{Float64,Vector{Dict{String,Float64}}}()
    for epsilon in EPSILONS
        token=replace(@sprintf("%.3g",epsilon),"."=>"p")
        rows=dict_csv(joinpath(STUDY_ROOT,"N360","mode_epsilon_$(token).csv"))
        mode_rows[epsilon]=filter(r->isfinite(r["Y"]),rows)
    end
    ygrid=collect(range(0,10;length=1001))
    fields=[("vG_scaled","v_G"),("vT_scaled","v_T"),
            ("vH_scaled","v_H / epsilon"),("vF_scaled","v_F / epsilon^2")]
    curves=Dict{Tuple{Float64,String},Vector{Float64}}()
    for epsilon in EPSILONS, (field,_) in fields
        rows=mode_rows[epsilon]
        x=[r["Y"] for r in rows]; y=[r[field] for r in rows]
        curves[(epsilon,field)]=interpolation(x,y,ygrid)
    end
    errors=NamedTuple[]
    reference=0.01
    for epsilon in (0.02,0.015), (field,_) in fields
        a=curves[(epsilon,field)]; b=curves[(reference,field)]
        absolute=norm(a-b)/sqrt(length(a))
        relative=norm(a-b)/max(norm(b),1e-14)
        push!(errors,(;epsilon,reference,field,absolute,relative))
    end
    open(joinpath(STUDY_ROOT,"mode_collapse_errors.csv"),"w") do io
        println(io,"epsilon,reference_epsilon,field,rms_difference,relative_l2_difference")
        for r in errors
            @printf(io,"%.12e,%.12e,%s,%.12e,%.12e\n",
                r.epsilon,r.reference,r.field,r.absolute,r.relative)
        end
    end
    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8,legendfontsize=8)
    panels=Any[]
    colors=Dict(0.02=>:navy,0.015=>:darkorange,0.01=>:seagreen)
    for (field,label) in fields
        panel=plot(;xlabel="Y = epsilon eta",ylabel=label,xlim=(0,10))
        for epsilon in EPSILONS
            plot!(panel,ygrid,curves[(epsilon,field)];color=colors[epsilon],marker=:none,
                  label=@sprintf("epsilon=%.3f",epsilon))
        end
        push!(panels,panel)
    end
    savefig(plot(panels...;layout=(2,2),size=(980,720)),
            joinpath(STUDY_ROOT,"fold_mode_outer_collapse_N360.png"))
end

function cusp_scaling()
    source=joinpath(ROOT,"work","results","consistent_fold_tail_continuation",
                    "apparent_cusp_migration.csv")
    rows=dict_csv(source)
    N=[r["N"] for r in rows]
    epsilon=[-r["Hinf"] for r in rows]
    coefficient=linear_fit(log.(N),log.(epsilon))
    exponent=-coefficient[2]
    prefactor=exp(coefficient[1])
    n2constant=mean(epsilon.*N.^2)
    open(joinpath(STUDY_ROOT,"cusp_resolution_scaling.csv"),"w") do io
        println(io,"N,epsilon_turn,N2_epsilon_turn,pairwise_power")
        for i in eachindex(N)
            pairpower=i==1 ? NaN : log(epsilon[i-1]/epsilon[i])/log(N[i]/N[i-1])
            @printf(io,"%.0f,%.12e,%.12e,%.12e\n",N[i],epsilon[i],N[i]^2*epsilon[i],pairpower)
        end
    end
    open(joinpath(STUDY_ROOT,"cusp_resolution_fit.txt"),"w") do io
        @printf(io,"loglog_power=%.12f\n",exponent)
        @printf(io,"loglog_prefactor=%.12f\n",prefactor)
        @printf(io,"mean_N2_epsilon=%.12f\n",n2constant)
        println(io,"scope=three existing apparent-turn points; numerical resolution boundary, not a physical cusp law")
    end
    nplot=collect(range(minimum(N),220;length=300))
    panel=scatter(N,epsilon;xscale=:log10,yscale=:log10,marker=:circle,
        label="apparent turn",xlabel="N",ylabel="epsilon_min = -Hinf_turn")
    plot!(panel,nplot,prefactor.*nplot.^(-exponent);marker=:none,label=@sprintf("fit N^{-%.2f}",exponent))
    plot!(panel,nplot,n2constant.*nplot.^(-2);marker=:none,linestyle=:dash,label="N^{-2}")
    savefig(panel,joinpath(STUDY_ROOT,"cusp_resolution_boundary.png"))
end

function main()
    mkpath(STUDY_ROOT)
    rows=tangent_rows()
    open(joinpath(STUDY_ROOT,"tangent_convergence.csv"),"w") do io
        println(io,"N,epsilon,dRo_depsilon,dTw_depsilon,dTw_dRo,tangent_residual")
        for r in sort(rows,by=q->(q.epsilon,q.N),rev=true)
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                r.N,r.epsilon,r.dRo,r.dTw,r.ratio,r.residual)
        end
    end
    fits=NamedTuple[]
    for N in DEGREES
        selected=sort(filter(r->r.N==N,rows),by=r->r.epsilon)
        epsilon=getfield.(selected,:epsilon)
        rofit=linear_fit(epsilon,getfield.(selected,:dRo))
        twfit=linear_fit(epsilon,getfield.(selected,:dTw))
        push!(fits,(;N,AR=rofit[1],AT=twfit[1],Ro_slope=rofit[2],Tw_slope=twfit[2]))
    end
    open(joinpath(STUDY_ROOT,"three_point_tangent_fits.csv"),"w") do io
        println(io,"N,A_R_linear,A_T_linear,dRo_tangent_slope,dTw_tangent_slope")
        for r in fits
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e\n",r.N,r.AR,r.AT,r.Ro_slope,r.Tw_slope)
        end
    end
    resolved_fits=NamedTuple[]
    for N in DEGREES
        selected=sort(filter(r->r.N==N && r.epsilon>=0.015,rows),by=r->r.epsilon)
        epsilon=getfield.(selected,:epsilon)
        rofit=linear_fit(epsilon,getfield.(selected,:dRo))
        twfit=linear_fit(epsilon,getfield.(selected,:dTw))
        push!(resolved_fits,(;N,AR=rofit[1],AT=twfit[1]))
    end
    open(joinpath(STUDY_ROOT,"resolved_window_tangent_fits.csv"),"w") do io
        println(io,"N,A_R_from_epsilon_0p02_0p015,A_T_from_epsilon_0p02_0p015")
        for r in resolved_fits
            @printf(io,"%d,%.12e,%.12e\n",r.N,r.AR,r.AT)
        end
    end
    open(joinpath(STUDY_ROOT,"resolution_spread.csv"),"w") do io
        println(io,"epsilon,field,N_min,N_max,minimum,maximum,relative_spread")
        for epsilon in EPSILONS, field in (:dRo,:dTw)
            selected=filter(r->r.epsilon==epsilon && r.N>=240,rows)
            values=getfield.(selected,field)
            relative=(maximum(values)-minimum(values))/abs(mean(values))
            @printf(io,"%.12e,%s,240,360,%.12e,%.12e,%.12e\n",
                epsilon,string(field),minimum(values),maximum(values),relative)
        end
    end
    default(;framestyle=:box,gridalpha=0.25,linewidth=2,marker=:circle,
            guidefontsize=10,tickfontsize=8,legendfontsize=8)
    panels=Any[]
    for (field,label,target) in ((:dRo,"dRo/d epsilon",-0.420342495108),
                                 (:dTw,"dTw/d epsilon",0.510660114111))
        panel=plot(;xlabel="epsilon=-Hinf",ylabel=label)
        for N in DEGREES
            selected=sort(filter(r->r.N==N,rows),by=r->r.epsilon)
            plot!(panel,getfield.(selected,:epsilon),getfield.(selected,field);
                  label="N=$N")
        end
        hline!(panel,[target];color=:black,linestyle=:dash,label="branch-fit limit")
        push!(panels,panel)
    end
    savefig(plot(panels...;layout=(1,2),size=(980,390)),
            joinpath(STUDY_ROOT,"tangent_convergence.png"))
    mode_collapse()
    cusp_scaling()
    println("Outputs: ",STUDY_ROOT)
end

main()

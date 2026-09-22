#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using DelimitedFiles
using Printf

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DEGREE=360
const EPSILONS=[0.02,0.015,0.01]
const PR=0.72; const GAMMA=1.0

function linear_interpolate(x,y,xnew,farvalue)
    finite=isfinite.(x)
    xf=x[finite]; yf=y[finite]
    output=similar(xnew)
    for i in eachindex(xnew)
        q=xnew[i]
        if !isfinite(q)
            output[i]=farvalue
        elseif q<=xf[1]
            output[i]=yf[1]
        elseif q>=xf[end]
            output[i]=farvalue
        else
            j=searchsortedlast(xf,q)
            weight=(q-xf[j])/(xf[j+1]-xf[j])
            output[i]=(1-weight)*yf[j]+weight*yf[j+1]
        end
    end
    output
end

function load_N320_source()
    directory=joinpath(ROOT,"work","results","consistent_epsilon_fold_chain","convergence_N320")
    table=readdlm(joinpath(directory,"epsilon_fold_table.csv"),',';skipstart=1)
    row=table[findfirst(==(0.02),table[:,1]),:]
    data=readdlm(joinpath(directory,"profiles","profile_epsilon_0p02.csv"),',';skipstart=1)
    row,data
end

function initial_N360_fold()
    row,data=load_N320_source()
    ro,tw,hinf=row[2],row[3],row[4]
    op=BEKIsothermal._operators(DEGREE,2.0,0.6,0.5)
    far=(-0.02,0.0,1.0,1.0)
    fields=zeros(4,length(op.x))
    for j in 1:4
        fields[j,:].=linear_interpolate(data[:,1],data[:,j+1],op.eta,far[j])
    end
    seed=ConsistentSolution(ro,2-ro-ro^2,GAMMA,PR,op,fields,tw,hinf,Inf,0)
    println("N=360: correcting interpolated base flow at epsilon=0.02 ...")
    base=solve_consistent_fixed_h_at_ro(-0.02,ro,seed;tolerance=1e-9)
    println("N=360: correcting fixed-Ro fold ...")
    fold=solve_consistent_fold(base;tolerance=1e-9)
    println("N=360: releasing Ro at fixed epsilon=0.02 ...")
    solve_consistent_fold_at_h(-0.02,fold;tolerance=1e-9)
end

function left_slope(fold)
    sol=fold.solution; op=sol.operators; state=BEKConsistent.pack(sol.fields)
    _,J=BEKConsistent._residual_jacobian(state,sol.Ro,sol.Tw,GAMMA,PR,op)
    rownorm=max.([norm(view(J,i,:)) for i in axes(J,1)],1e-14)
    decomp=svd(J./rownorm)
    left=(1 ./rownorm).*decomp.U[:,end]; left./=norm(left)
    rro=BEKConsistent._residual_ro_derivative(state,sol.Ro,sol.Tw,GAMMA,PR,op)
    rtw=BEKConsistent._tw_derivative(state,op)
    -dot(left,rro)/dot(left,rtw)
end

function save_point(output,epsilon,fold,tangent)
    sol=fold.solution; n=length(sol.operators.x)
    mode=BEKConsistent.unpack(fold.nullvector,n)
    maxg=maximum(abs,view(mode,3,:))
    open(joinpath(output,"mode_epsilon_$(replace(@sprintf("%.3g",epsilon),"."=>"p")).csv"),"w") do io
        println(io,"eta,Y,vH_scaled,vF_scaled,vG_scaled,vT_scaled")
        for i in 1:n
            @printf(io,"%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                sol.operators.eta[i],epsilon*sol.operators.eta[i],
                mode[1,i]/(maxg*epsilon),mode[2,i]/(maxg*epsilon^2),
                mode[3,i]/maxg,mode[4,i]/maxg)
        end
    end
end

function main()
    output=joinpath(ROOT,"work","results","tangent_convergence_mode_collapse","N360")
    mkpath(output)
    fold=initial_N360_fold()
    rows=NamedTuple[]
    for epsilon in EPSILONS
        if abs(fold.solution.Hinf+epsilon)>1e-12
            println("N=360: advancing fold to epsilon=$epsilon ...")
            fold=solve_consistent_fold_at_h(-epsilon,fold;tolerance=1e-9)
        end
        tangent=fold_hinf_tangent(fold)
        slope=left_slope(fold)
        push!(rows,(;epsilon,Ro=fold.solution.Ro,Tw=fold.solution.Tw,
            fold_residual=fold.residual,dRo=tangent.Ro,dTw=tangent.Tw,
            ratio=tangent.Tw/tangent.Ro,left_slope=slope,
            tangent_residual=tangent.residual))
        save_point(output,epsilon,fold,tangent)
        @printf("epsilon=%.3f Ro=%+.10f dRo/deps=%+.9f dTw/deps=%+.9f ratio=%+.9f\n",
            epsilon,fold.solution.Ro,tangent.Ro,tangent.Tw,tangent.Tw/tangent.Ro)
    end
    open(joinpath(output,"tangent_table.csv"),"w") do io
        println(io,"N,epsilon,Ro,Tw,fold_residual,dRo_depsilon,dTw_depsilon,dTw_dRo,left_slope,tangent_residual")
        for r in rows
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                DEGREE,r.epsilon,r.Ro,r.Tw,r.fold_residual,r.dRo,r.dTw,
                r.ratio,r.left_slope,r.tangent_residual)
        end
    end
end

main()

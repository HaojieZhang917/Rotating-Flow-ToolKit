#!/usr/bin/env julia

include(joinpath(@__DIR__,"..","src","BEKIsothermal.jl"))
include(joinpath(@__DIR__,"..","src","BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra
using DelimitedFiles
using Printf

function option(prefix,default)
    for argument in ARGS
        startswith(argument,prefix) && return split(argument,"=";limit=2)[2]
    end
    default
end

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DEGREE=parse(Int,option("--degree=","240"))
const SOURCE_EPS=0.02
const TARGET_EPS=0.015
const PR=0.72; const GAMMA=1.0
const RUN_NAMES=Dict(200=>"convergence_N200",240=>"convergence_N240",
                     280=>"convergence_N280",320=>"convergence_N320")

function numeric_csv(path)
    lines=filter(!isempty,strip.(readlines(path)))
    header=Symbol.(split(lines[1],','))
    [NamedTuple{Tuple(header)}(Tuple(parse.(Float64,split(line,',')))) for line in lines[2:end]]
end

token(epsilon)=replace(@sprintf("%.4g",epsilon),"."=>"p","-"=>"m")

function load_source_fold()
    run=RUN_NAMES[DEGREE]
    directory=joinpath(ROOT,"work","results","consistent_epsilon_fold_chain",run)
    row=only(filter(r->r.epsilon==SOURCE_EPS,
        numeric_csv(joinpath(directory,"epsilon_fold_table.csv"))))
    data=readdlm(joinpath(directory,"profiles","profile_epsilon_$(token(SOURCE_EPS)).csv"),',';skipstart=1)
    op=BEKIsothermal._operators(DEGREE,2.0,0.6,0.5)
    fields=permutedims(data[:,2:5])
    state=BEKConsistent.pack(fields)
    residual,jacobian=BEKConsistent._residual_jacobian(state,row.Ro_f,row.Tw_f,GAMMA,PR,op)
    solution=ConsistentSolution(row.Ro_f,2-row.Ro_f-row.Ro_f^2,GAMMA,PR,op,
        fields,row.Tw_f,row.Hinf,norm(residual,Inf),0)
    rownorm=max.([norm(view(jacobian,i,:)) for i in axes(jacobian,1)],1e-14)
    decomposition=svd(jacobian./rownorm)
    vector=decomposition.V[:,end]
    vector./=norm(vector)
    vector[length(op.x)]>0 && (vector.*=-1)
    ConsistentFold(solution,vector,norm(residual,Inf),0)
end

function left_mode_and_slope(fold)
    solution=fold.solution; op=solution.operators
    state=BEKConsistent.pack(solution.fields)
    _,jacobian=BEKConsistent._residual_jacobian(state,solution.Ro,
        solution.Tw,GAMMA,PR,op)
    rownorm=max.([norm(view(jacobian,i,:)) for i in axes(jacobian,1)],1e-14)
    decomposition=svd(jacobian./rownorm)
    left=(1 ./rownorm).*decomposition.U[:,end]
    left./=norm(left)
    rro=BEKConsistent._residual_ro_derivative(state,solution.Ro,
        solution.Tw,GAMMA,PR,op)
    rtw=BEKConsistent._tw_derivative(state,op)
    slope=-dot(left,rro)/dot(left,rtw)
    slope,decomposition.S[end],decomposition.S[end]/decomposition.S[1]
end

function main()
    output=joinpath(ROOT,"work","results","tangent_convergence_mode_collapse",
                    "N$(DEGREE)")
    mkpath(output)
    println("N=$DEGREE: correcting existing epsilon=0.02 fold to epsilon=0.015 ...")
    source=load_source_fold()
    fold=solve_consistent_fold_at_h(-TARGET_EPS,source;tolerance=1e-9)
    tangent=fold_hinf_tangent(fold)
    slope,sigma_min,sigma_ratio=left_mode_and_slope(fold)
    solution=fold.solution; n=length(solution.operators.x)
    mode=BEKConsistent.unpack(fold.nullvector,n)
    maxg=maximum(abs,view(mode,3,:))
    open(joinpath(output,"fold_and_tangent.csv"),"w") do io
        println(io,"N,epsilon,Ro,Tw,Hinf,fold_residual,dRo_depsilon,dTw_depsilon,dTw_dRo,left_slope,dHinf_depsilon,tangent_residual,scaled_sigma_min,scaled_sigma_ratio")
        @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
            DEGREE,TARGET_EPS,solution.Ro,solution.Tw,solution.Hinf,fold.residual,
            tangent.Ro,tangent.Tw,tangent.Tw/tangent.Ro,slope,tangent.dHinf,
            tangent.residual,sigma_min,sigma_ratio)
    end
    tangent_fields=BEKConsistent.unpack(tangent.field,n)
    open(joinpath(output,"mode_epsilon_0p015.csv"),"w") do io
        println(io,"eta,Y,vH,vF,vG,vT,vH_scaled,vF_scaled,vG_scaled,vT_scaled,H_epsilon,F_epsilon,G_epsilon,T_epsilon")
        for i in 1:n
            eta=solution.operators.eta[i]; Y=TARGET_EPS*eta
            @printf(io,"%.12e,%.12e",eta,Y)
            for value in mode[:,i]; @printf(io,",%.12e",value); end
            scaled=(mode[1,i]/(maxg*TARGET_EPS),
                    mode[2,i]/(maxg*TARGET_EPS^2),
                    mode[3,i]/maxg,mode[4,i]/maxg)
            for value in scaled; @printf(io,",%.12e",value); end
            for value in tangent_fields[:,i]; @printf(io,",%.12e",value); end
            println(io)
        end
    end
    @printf("Ro=%+.12f Tw=%.12f dRo/deps=%+.10f dTw/deps=%+.10f ratio=%+.10f residual=%.2e\n",
        solution.Ro,solution.Tw,tangent.Ro,tangent.Tw,tangent.Tw/tangent.Ro,
        tangent.residual)
end

main()

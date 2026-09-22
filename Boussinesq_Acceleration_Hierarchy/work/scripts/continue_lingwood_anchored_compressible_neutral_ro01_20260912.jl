using BSplineKit
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using PyCall

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const TOOLKIT_ROOT = normpath(joinpath(ROOT,".."))
const BACKEND = joinpath(TOOLKIT_ROOT,"RotatingDiskFlow","src","CRD_STA.jl")
const FRONTEND = joinpath(ROOT,"work","src","CompressibleBEKLinearStability.jl")
const OUTDIR = joinpath(ROOT,"work","results",
    "compressible_bek_lingwood_anchored_neutral_v3_20260915")

function load_solver()
    backend_text = replace(read(BACKEND,String),
        "using DifferentialEquations"=>"",
        "using BoundaryValueDiffEq"=>"")
    include_string(Main,backend_text,BACKEND*"#anchored-neutral")
    frontend_text = replace(read(FRONTEND,String),
        "using DifferentialEquations"=>"",
        "using BoundaryValueDiffEq"=>"")
    old_block = """const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at \$BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new_block = """const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CRD = Main"""
    occursin(old_block,frontend_text) || error("frontend loader block changed")
    include_string(Main,replace(frontend_text,old_block=>new_block),FRONTEND*"#anchored-neutral")
    Main.CompressibleBEKLinearStability
end

const CBEK = load_solver()
const MACH_PATH = [0.01,0.05,0.10,0.15,0.20,0.25,0.30]
const FLOWS = [
    (flow="Ekman",Ro=0.0,R0=116.3,beta0=0.137,alpha0=0.528,beta_target=0.132,
     beta_path=[0.136,0.135,0.134,0.133,0.132],
     lower=[110.0,100.0,90.0,80.0,70.0],upper=[125.0,140.0,160.0,180.0]),
    (flow="Bodewadt",Ro=1.0,R0=27.4,beta0=0.115,alpha0=0.487,beta_target=0.127,
     beta_path=[0.117,0.119,0.121,0.123,0.125,0.127],
     lower=[25.0,22.5,20.0,17.5,15.0],upper=[30.0,35.0,40.0,50.0,60.0]),
]

function make_case(flow;Mr,R,beta,target,N=49)
    p=CBEK.PhysicalParameters(Mr=Mr,Tw=1.0,Ro=flow.Ro,R=R,beta=beta,
        omega=0.0,gamma=1.4,Pr=0.72)
    n=CBEK.NumericalParameters(N=N,domain=:rational,target=ComplexF64(target),
        neigs=1,tol=1e-9,maxit=600,balance=true,balance_iterations=2,
        regularization=0.0)
    CBEK.prepare_case(p,n;allow_unverified_ro=true)
end

function write_rows(path,rows)
    isempty(rows) && return
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(row->println(io,join(values(row),',')),rows)
    end
end

function solve_point(flow;Mr,R,beta,target,vector=nothing,N=49,stage,index)
    case=make_case(flow;Mr=Mr,R=R,beta=beta,target=target,N=N)
    result=CBEK.solve_mode(case;initial_vector=vector)
    alpha=first(result.alpha); overlap=first(result.overlaps)
    row=(flow=flow.flow,stage=stage,index=index,Ro=flow.Ro,Mr=Mr,N=N,R=R,
        beta=beta,alpha_real=real(alpha),alpha_imag=imag(alpha),overlap=overlap,
        backward_error=first(result.backward_errors),target_real=real(target),
        target_imag=imag(target),evidence_scope=String(case.evidence_scope))
    @printf("%-8s %-10s %2d N=%d beta=%.3f R=%8.4f alpha=%.10f%+.10fi ov=%s\n",
        flow.flow,stage,index,N,beta,R,real(alpha),imag(alpha),
        isnan(overlap) ? "seed" : @sprintf("%.6f",overlap))
    (row=row,alpha=alpha,vector=copy(view(result.vectors,:,1)))
end

function append_point!(rows,point)
    push!(rows,point.row)
    write_rows(joinpath(OUTDIR,"continuation_points.csv"),rows)
    point
end

function advance(flow,start,values,parameter,rows,stage)
    point=start
    for (index,value) in enumerate(values)
        kwargs=parameter==:Mr ? (Mr=value,R=flow.R0,beta=flow.beta0) :
            parameter==:beta ? (Mr=0.3,R=flow.R0,beta=value) :
            (Mr=0.3,R=value,beta=flow.beta_target)
        point=solve_point(flow;kwargs...,target=point.alpha,vector=point.vector,
            stage=stage,index=index)
        append_point!(rows,point)
    end
    point
end

function radial_points(flow,seed,values,rows,stage)
    points=[seed]; point=seed
    for (index,R) in enumerate(values)
        point=solve_point(flow;Mr=0.3,R=R,beta=flow.beta_target,
            target=point.alpha,vector=point.vector,stage=stage,index=index)
        append_point!(rows,point); push!(points,point)
    end
    points
end

function brackets(points)
    [(points[j],points[j+1]) for j in 1:length(points)-1 if
        signbit(imag(points[j].alpha)) != signbit(imag(points[j+1].alpha))]
end

function refine!(flow,left,right,rows,candidates)
    left.row.R>right.row.R && ((left,right)=(right,left))
    point=left
    for index in 1:20
        fl,fr=imag(left.alpha),imag(right.alpha)
        R=(left.row.R*fr-right.row.R*fl)/(fr-fl)
        width=right.row.R-left.row.R
        (!isfinite(R) || R<=left.row.R+0.1width || R>=right.row.R-0.1width) &&
            (R=(left.row.R+right.row.R)/2)
        seed=abs(R-left.row.R)<=abs(R-right.row.R) ? left : right
        point=solve_point(flow;Mr=0.3,R=R,beta=flow.beta_target,
            target=seed.alpha,vector=seed.vector,stage="refine",index=index)
        append_point!(rows,point)
        abs(imag(point.alpha))<1e-8 && break
        signbit(imag(left.alpha)) != signbit(imag(point.alpha)) ?
            (right=point) : (left=point)
    end
    for N in (49,59,69)
        check=N==49 ? point : solve_point(flow;Mr=0.3,R=point.row.R,
            beta=flow.beta_target,target=point.alpha,N=N,stage="grid",index=N)
        push!(candidates,(flow=flow.flow,Ro=flow.Ro,Mr=0.3,N=N,R=point.row.R,
            beta=flow.beta_target,alpha_real=real(check.alpha),
            alpha_imag=imag(check.alpha),backward_error=check.row.backward_error,
            evidence_scope=check.row.evidence_scope))
        write_rows(joinpath(OUTDIR,"neutral_candidates.csv"),candidates)
    end
end

mkpath(OUTDIR)
rows=NamedTuple[]; candidates=NamedTuple[]; summaries=NamedTuple[]
for flow in FLOWS
    seed=solve_point(flow;Mr=0.01,R=flow.R0,beta=flow.beta0,
        target=ComplexF64(flow.alpha0),stage="seed",index=1)
    append_point!(rows,seed)
    mach=advance(flow,seed,MACH_PATH[2:end],:Mr,rows,"mach")
    beta=advance(flow,mach,flow.beta_path,:beta,rows,"beta")
    down=radial_points(flow,beta,flow.lower,rows,"R_down")
    up=radial_points(flow,beta,flow.upper,rows,"R_up")
    found=vcat(brackets(down),brackets(up))
    push!(summaries,(flow=flow.flow,Ro=flow.Ro,Mr=0.3,N=49,
        beta=flow.beta_target,sign_change_count=length(found),
        anchored_alpha_real=real(beta.alpha),anchored_alpha_imag=imag(beta.alpha)))
    write_rows(joinpath(OUTDIR,"scan_summary.csv"),summaries)
    foreach(pair->refine!(flow,pair[1],pair[2],rows,candidates),found)
end

open(joinpath(OUTDIR,"metadata.txt"),"w") do io
    println(io,"study=Lingwood_anchored_corrected_compressible_fixed_beta_neutral")
    println(io,"path=low_Mach_Lingwood_mode_to_Mr_0.3_to_supplied_beta_to_neutral_R")
    println(io,"operator_source=$(FRONTEND)")
    println(io,"legacy_backend_preserved=$(BACKEND)")
end
println("Wrote $(OUTDIR)")

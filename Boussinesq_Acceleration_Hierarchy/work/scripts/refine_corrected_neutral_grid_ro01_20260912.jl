using BoundaryValueDiffEq
using BSplineKit
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using PyCall

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const TOOLKIT_ROOT=normpath(joinpath(ROOT,".."))
const BACKEND=joinpath(TOOLKIT_ROOT,"RotatingDiskFlow","src","CRD_STA.jl")
const FRONTEND=joinpath(ROOT,"work","src","CompressibleBEKLinearStability.jl")
const OUTDIR=joinpath(ROOT,"work","results",
    "compressible_bek_lingwood_anchored_neutral_20260912")

function load_solver()
    backend_text=replace(read(BACKEND,String),
        "using DifferentialEquations"=>"using BoundaryValueDiffEq")
    include_string(Main,backend_text,BACKEND*"#neutral-grid")
    frontend_text=replace(read(FRONTEND,String),
        "using DifferentialEquations"=>"using BoundaryValueDiffEq")
    old_block="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at \$BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new_block="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CRD = Main"""
    include_string(Main,replace(frontend_text,old_block=>new_block),FRONTEND*"#neutral-grid")
    Main.CompressibleBEKLinearStability
end
const CBEK=load_solver()

const CASES=[
    (flow="Ekman",Ro=0.0,beta=0.132,left=116.3,right=125.0,target=0.49804+0im),
    (flow="Bodewadt",Ro=1.0,beta=0.127,left=30.0,right=35.0,target=0.61692+0im),
]

function make_case(c,N,R,target)
    p=CBEK.PhysicalParameters(Mr=0.3,Tw=1.0,Ro=c.Ro,R=R,beta=c.beta,
        omega=0.0,gamma=1.4,Pr=0.72)
    n=CBEK.NumericalParameters(N=N,domain=:rational,target=ComplexF64(target),
        neigs=1,tol=1e-9,maxit=700,balance=true,balance_iterations=2,
        regularization=0.0)
    CBEK.prepare_case(p,n;allow_unverified_ro=true)
end

function solve_at(c,N,R,target,vector=nothing)
    result=CBEK.solve_mode(make_case(c,N,R,target);initial_vector=vector)
    (R=R,alpha=first(result.alpha),vector=copy(view(result.vectors,:,1)),
     overlap=first(result.overlaps),backward_error=first(result.backward_errors))
end

function refine(c,N)
    left=solve_at(c,N,c.left,c.target)
    right=solve_at(c,N,c.right,left.alpha,left.vector)
    imag(left.alpha)*imag(right.alpha)<=0 || error("no $(c.flow) bracket at N=$N")
    point=left
    for index in 1:22
        fl,fr=imag(left.alpha),imag(right.alpha)
        R=(left.R*fr-right.R*fl)/(fr-fl)
        width=right.R-left.R
        (!isfinite(R) || R<=left.R+0.1width || R>=right.R-0.1width) &&
            (R=(left.R+right.R)/2)
        seed=abs(R-left.R)<=abs(R-right.R) ? left : right
        point=solve_at(c,N,R,seed.alpha,seed.vector)
        @printf("%-8s N=%d iter=%2d R=%.9f alpha=%.11f%+.3ei ov=%.7f\n",
            c.flow,N,index,R,real(point.alpha),imag(point.alpha),point.overlap)
        abs(imag(point.alpha))<1e-9 && break
        signbit(imag(left.alpha)) != signbit(imag(point.alpha)) ?
            (right=point) : (left=point)
    end
    (flow=c.flow,Ro=c.Ro,Mr=0.3,N=N,R=point.R,beta=c.beta,
     alpha_real=real(point.alpha),alpha_imag=imag(point.alpha),
     backward_error=point.backward_error,final_overlap=point.overlap,
     evidence_scope="derived_bek_candidate")
end

rows=[refine(c,N) for c in CASES for N in (49,59,69)]
open(joinpath(OUTDIR,"neutral_grid_convergence.csv"),"w") do io
    println(io,join(keys(first(rows)),','))
    foreach(row->println(io,join(values(row),',')),rows)
end
println("Wrote $(joinpath(OUTDIR,"neutral_grid_convergence.csv"))")

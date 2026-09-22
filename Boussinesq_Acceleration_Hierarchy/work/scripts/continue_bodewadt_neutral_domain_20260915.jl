using BSplineKit
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using PyCall

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const TOOLKIT_ROOT=normpath(joinpath(ROOT,".."))
const BACKEND=joinpath(TOOLKIT_ROOT,"RotatingDiskFlow","src","CRD_STA.jl")
const FRONTEND=joinpath(ROOT,"work","src","CompressibleBEKLinearStability.jl")
const OUTDIR=joinpath(ROOT,"work","results","compressible_bek_bodewadt_neutral_domain_20260915")

function load_solver()
    b=replace(read(BACKEND,String),"using DifferentialEquations"=>"","using BoundaryValueDiffEq"=>"")
    include_string(Main,b,BACKEND*"#bode-domain")
    f=replace(read(FRONTEND,String),"using DifferentialEquations"=>"","using BoundaryValueDiffEq"=>"")
    old="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at \$BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CRD = Main"""
    occursin(old,f) || error("frontend loader block changed")
    include_string(Main,replace(f,old=>new),FRONTEND*"#bode-domain")
    Main.CompressibleBEKLinearStability
end
const CBEK=load_solver()

function solve_point(L,R,target,vector)
    p=CBEK.PhysicalParameters(Mr=0.3,Tw=1.0,Ro=1.0,R=R,beta=0.127,
        omega=0.0,gamma=1.4,Pr=0.72)
    n=CBEK.NumericalParameters(N=49,domain=:rational,ymax=L,
        target=ComplexF64(target),neigs=1,tol=1e-9,maxit=600)
    c=CBEK.prepare_case(p,n;allow_unverified_ro=true)
    s=CBEK.solve_mode(c;initial_vector=vector)
    (R=R,alpha=first(s.alpha),vector=copy(view(s.vectors,:,1)),backward=first(s.backward_errors))
end

mkpath(OUTDIR)
rows=NamedTuple[]
summary=NamedTuple[]
for L in (40.0,60.0,80.0)
    points=NamedTuple[]
    target=0.50+0im; vector=nothing
    for R in (35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,80.0)
        p=solve_point(L,R,target,vector); push!(points,p)
        @printf("L=%4.0f R=%6.1f alpha_i=%+.6e\\n",L,R,imag(p.alpha))
        push!(rows,(L=L,R=R,alpha_real=real(p.alpha),alpha_imag=imag(p.alpha),
            backward_error=p.backward,stage="scan"))
        target=p.alpha; vector=p.vector
    end
    brackets=Tuple{Any,Any}[]
    for i in 1:length(points)-1
        signbit(imag(points[i].alpha)) != signbit(imag(points[i+1].alpha)) &&
            push!(brackets,(points[i],points[i+1]))
    end
    if isempty(brackets)
        push!(summary,(L=L,R_neutral=NaN,alpha_real=NaN,alpha_imag=NaN,
            bracket_left=NaN,bracket_right=NaN))
        continue
    end
    left,right=first(brackets)
    for k in 1:14
        fl,fr=imag(left.alpha),imag(right.alpha)
        R=(left.R*fr-right.R*fl)/(fr-fl)
        R<=left.R+0.1*(right.R-left.R) && (R=(left.R+right.R)/2)
        R>=right.R-0.1*(right.R-left.R) && (R=(left.R+right.R)/2)
        seed=abs(R-left.R)<abs(R-right.R) ? left : right
        p=solve_point(L,R,seed.alpha,seed.vector)
        push!(rows,(L=L,R=R,alpha_real=real(p.alpha),alpha_imag=imag(p.alpha),
            backward_error=p.backward,stage="refine"))
        if abs(imag(p.alpha))<1e-9
            left=right=p; break
        elseif signbit(imag(left.alpha)) != signbit(imag(p.alpha))
            right=p
        else
            left=p
        end
    end
    root=abs(imag(left.alpha))<abs(imag(right.alpha)) ? left : right
    for N in (49,59,69)
        p=CBEK.PhysicalParameters(Mr=0.3,Tw=1.0,Ro=1.0,R=root.R,beta=0.127,
            omega=0.0,gamma=1.4,Pr=0.72)
        n=CBEK.NumericalParameters(N=N,domain=:rational,ymax=L,
            target=root.alpha,neigs=1,tol=1e-9,maxit=600)
        c=CBEK.prepare_case(p,n;allow_unverified_ro=true)
        s=CBEK.solve_mode(c)
        push!(rows,(L=L,R=root.R,alpha_real=real(first(s.alpha)),
            alpha_imag=imag(first(s.alpha)),backward_error=first(s.backward_errors),
            stage="grid_$N"))
    end
    push!(summary,(L=L,R_neutral=root.R,alpha_real=real(root.alpha),
        alpha_imag=imag(root.alpha),bracket_left=first(brackets)[1].R,
        bracket_right=first(brackets)[2].R))
    @printf("L=%4.0f neutral R=%.8f alpha=%.8f%+.3ei\n",L,root.R,
        real(root.alpha),imag(root.alpha))
end
function write_rows(path,rs)
    open(path,"w") do io
        println(io,join(keys(first(rs)),',')); foreach(r->println(io,join(values(r),',')),rs)
    end
end
write_rows(joinpath(OUTDIR,"continuation.csv"),rows)
write_rows(joinpath(OUTDIR,"summary.csv"),summary)
open(joinpath(OUTDIR,"metadata.txt"),"w") do io
    println(io,"study=Bodewadt_Mr03_neutral_domain_convergence")
    println(io,"domains=40,60,80; scan_N=49; grid_checks=49,59,69")
    println(io,"fixed_beta=0.127; Tw=1; Mr=0.3; Ro=1")
end
println("Wrote $(OUTDIR)")

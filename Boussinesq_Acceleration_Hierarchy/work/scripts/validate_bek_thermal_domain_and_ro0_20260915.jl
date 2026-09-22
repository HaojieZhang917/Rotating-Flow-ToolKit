using BSplineKit
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using PyCall

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const TOOLKIT_ROOT=normpath(joinpath(ROOT,".."))
const BACKEND=joinpath(TOOLKIT_ROOT,"RotatingDiskFlow","src","CRD_STA.jl")
const FRONTEND=joinpath(ROOT,"work","src","CompressibleBEKLinearStability.jl")
const RESIDUAL_SOURCE=joinpath(ROOT,"work","src","CompressibleBEKBaseflowResiduals.jl")
const OUTDIR=joinpath(ROOT,"work","results","compressible_bek_thermal_domain_ro0_20260915")

function load_solver()
    backend_text=replace(read(BACKEND,String),
        "using DifferentialEquations"=>"", "using BoundaryValueDiffEq"=>"")
    include_string(Main,backend_text,BACKEND*"#thermal-domain")
    frontend_text=replace(read(FRONTEND,String),
        "using DifferentialEquations"=>"", "using BoundaryValueDiffEq"=>"")
    old_block="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at \$BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new_block="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CRD = Main"""
    occursin(old_block,frontend_text) || error("frontend loader block changed")
    include_string(Main,replace(frontend_text,old_block=>new_block),FRONTEND*"#thermal-domain")
    Main.CompressibleBEKLinearStability
end

const CBEK=load_solver()
include(RESIDUAL_SOURCE)
const BRES=CompressibleBEKBaseflowResiduals

function sampled(source,y)
    fs=BSplineKit.interpolate(source.y,source.f,BSplineOrder(4))
    qs=BSplineKit.interpolate(source.y,source.q,BSplineOrder(4))
    f=[fs(x) for x in y]; q=[qs(x) for x in y]
    f[1]=0.0; f[end]=0.0; q[1]=1.0; q[end]=0.0
    f,q
end
function value(source,point,field)
    BSplineKit.interpolate(source.y,getfield(source,field),BSplineOrder(4))(point)
end
function write_rows(path,rows)
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(r->println(io,join(values(r),',')),rows)
    end
end

mkpath(OUTDIR)
rows=NamedTuple[]
for flow in ((name="Bodewadt",Ro=1.0,R=27.4,beta=0.115,target=0.487),
             (name="Ekman",Ro=0.0,R=116.3,beta=0.137,target=0.528)), L in (40.0,60.0,80.0)
    p=CBEK.PhysicalParameters(Mr=0.3,Tw=1.0,Ro=flow.Ro,R=flow.R,beta=0.0,
        omega=0.0,gamma=1.4,Pr=0.72)
    n=CBEK.NumericalParameters(N=80,domain=:finite,ymax=L,neigs=1)
    case=CBEK.prepare_case(p,n;allow_unverified_ro=true)
    source=CBEK.source_profiles(flow.Ro,p.Pr;ymax=L)
    f,q=sampled(source,case.y)
    check=BRES.baseflow_primitive_residuals(case.y,case.D,case.D2,case.F,case.G,
        case.H./case.T,case.T,case.rho,f,q;Ro=flow.Ro,Mr=p.Mr,Tw=p.Tw,
        gamma=p.gamma,Pr=p.Pr,max_y=L-1.0)
    points=(1.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0)
    qvals=[x<L ? value(source,x,:q) : NaN for x in points]
    mode_n=CBEK.NumericalParameters(N=49,domain=:rational,ymax=L,
        target=ComplexF64(flow.target),neigs=1,tol=1e-9)
    mode_p=CBEK.PhysicalParameters(Mr=0.3,Tw=1.0,Ro=flow.Ro,R=flow.R,
        beta=flow.beta,omega=0.0,gamma=1.4,Pr=0.72)
    mode_case=CBEK.prepare_case(mode_p,mode_n;allow_unverified_ro=true)
    mode_result=CBEK.solve_mode(mode_case)
    push!(rows,(flow=flow.name,Ro=flow.Ro,L=L,N=80,
        primitive_energy_abs=check.metrics.primitive_energy.absolute_linf,
        primitive_energy_term_ratio=check.metrics.primitive_energy.term_relative_linf,
        q_y1=qvals[1],q_y10=qvals[2],q_y20=qvals[3],q_y30=qvals[4],
        q_y40=qvals[5],q_y50=qvals[6],q_y60=qvals[7],q_y70=qvals[8],
        q_endminus1=value(source,L-1.0,:q),T_max=maximum(case.T),
        mode_alpha_real=real(first(mode_result.alpha)),
        mode_alpha_imag=imag(first(mode_result.alpha)),
        mode_backward_error=first(mode_result.backward_errors)))
end
write_rows(joinpath(OUTDIR,"domain_sweep.csv"),rows)
open(joinpath(OUTDIR,"metadata.txt"),"w") do io
    println(io,"study=compressible_BEK_thermal_domain_and_Ro0_degeneracy")
    println(io,"date=2026-09-15")
    println(io,"domains=40,60,80; finite Chebyshev N=80")
    println(io,"flows=Bodewadt Ro=1 and Ekman Ro=0")
end
println("Wrote $(OUTDIR)")

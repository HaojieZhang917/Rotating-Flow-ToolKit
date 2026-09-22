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

function load_solver()
    backend_text=replace(read(BACKEND,String),
        "using DifferentialEquations"=>"using BoundaryValueDiffEq")
    include_string(Main,backend_text,BACKEND*"#lingwood-gate")
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
    include_string(Main,replace(frontend_text,old_block=>new_block),FRONTEND*"#lingwood-gate")
    Main.CompressibleBEKLinearStability
end
const CBEK=load_solver()

function make_case(ref;Mr,N,target,neigs=1)
    p=CBEK.PhysicalParameters(Mr=Mr,Tw=1.0,Ro=ref.Ro,R=ref.R,
        beta=ref.beta,omega=0.0,gamma=1.4,Pr=0.72)
    n=CBEK.NumericalParameters(N=N,domain=:rational,target=target,neigs=neigs,
        tol=1e-9,maxit=600,balance=true,balance_iterations=2,regularization=0.0)
    CBEK.prepare_case(p,n;allow_unverified_ro=true)
end

function write_rows(path,rows)
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(row->println(io,join(values(row),',')),rows)
    end
end

const LINGWOOD_OUTDIR = joinpath(ROOT,"work","results",
    "compressible_bek_corrected_lingwood_gate_v2_20260912")
const LINGWOOD_REFERENCES = [
    (flow="von_Karman",Ro=-1.0,R=290.1,beta=0.077,alpha=0.381),
    (flow="Ekman",Ro=0.0,R=116.3,beta=0.137,alpha=0.528),
    (flow="Bodewadt",Ro=1.0,R=27.4,beta=0.115,alpha=0.487),
]

mkpath(LINGWOOD_OUTDIR)
lingwood_rows = NamedTuple[]
for ref in LINGWOOD_REFERENCES, N in (49,59)
    case = make_case(ref;Mr=0.01,N=N,target=ComplexF64(ref.alpha),neigs=1)
    result = CBEK.solve_mode(case)
    alpha = first(result.alpha)
    push!(lingwood_rows,(
        flow=ref.flow,Ro=ref.Ro,Mr=0.01,N=N,R=ref.R,beta=ref.beta,
        alpha_reference=ref.alpha,alpha_real=real(alpha),alpha_imag=imag(alpha),
        distance_to_reference=abs(alpha-ref.alpha),
        backward_error=first(result.backward_errors),
        evidence_scope=String(case.evidence_scope),
    ))
    @printf("Lingwood %-11s N=%d alpha=%.11f%+.11fi distance=%.3e\n",
        ref.flow,N,real(alpha),imag(alpha),abs(alpha-ref.alpha))
end
write_rows(joinpath(LINGWOOD_OUTDIR,"low_mach_lingwood_gate.csv"),lingwood_rows)
println("Wrote $(LINGWOOD_OUTDIR)")

using BSplineKit
using DelimitedFiles
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using PyCall

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TOOLKIT_ROOT = normpath(joinpath(ROOT, ".."))
const BACKEND = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "CRD_STA.jl")
const FRONTEND = joinpath(ROOT, "work", "src", "CompressibleBEKLinearStability.jl")
const OUTDIR = joinpath(ROOT, "work", "results",
    "compressible_bek_corrected_fixed_beta_neutral_v3_20260915")

function load_solver()
    backend_text = replace(read(BACKEND, String),
        "using DifferentialEquations" => "",
        "using BoundaryValueDiffEq" => "")
    include_string(Main, backend_text, BACKEND * "#minimal-import-neutral-continuation")
    frontend_text = replace(read(FRONTEND, String),
        "using DifferentialEquations" => "",
        "using BoundaryValueDiffEq" => "")
    old_block = """const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at \$BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new_block = """const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CRD = Main"""
    occursin(old_block, frontend_text) || error("frontend loader block changed")
    include_string(Main, replace(frontend_text, old_block => new_block),
        FRONTEND * "#minimal-import-neutral-continuation")
    return Main.CompressibleBEKLinearStability
end

const CBEK = load_solver()
const FLOWS = [
    (flow="Ekman", Ro=0.0, R=116.25, beta=0.132, alpha_ref=0.2775,
     lower=[110.0,100.0,90.0,80.0,70.0,60.0,50.0,40.0],
     upper=[125.0,140.0,160.0,180.0,210.0,250.0]),
    (flow="Bodewadt", Ro=1.0, R=27.26, beta=0.127, alpha_ref=0.5231,
     lower=[25.0,22.5,20.0,17.5,15.0,12.5,10.0],
     upper=[30.0,35.0,40.0,50.0,60.0,75.0]),
]
const MACH_PATH = [0.01,0.05,0.10,0.15,0.20,0.25,0.30]

function make_case(flow,Mr,R,target,N; neigs=1)
    p = CBEK.PhysicalParameters(
        Mr=Mr,Tw=1.0,Ro=flow.Ro,R=R,beta=flow.beta,
        omega=0.0,gamma=1.4,Pr=0.72,
    )
    n = CBEK.NumericalParameters(
        N=N,domain=:rational,target=ComplexF64(target),neigs=neigs,
        tol=1e-9,maxit=600,balance=true,balance_iterations=2,
        regularization=0.0,
    )
    return CBEK.prepare_case(p,n;allow_unverified_ro=true)
end

function write_rows(path,rows)
    isempty(rows) && return
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(row -> println(io,join(values(row),',')),rows)
    end
end

function solve_point(flow;Mr,R,N,target,previous_vector=nothing,stage,index)
    case = make_case(flow,Mr,R,target,N)
    result = CBEK.solve_mode(case;initial_vector=previous_vector)
    alpha = first(result.alpha)
    overlap = first(result.overlaps)
    row = (
        flow=flow.flow,stage=stage,index=index,Ro=flow.Ro,Mr=Mr,N=N,R=R,
        beta=flow.beta,alpha_real=real(alpha),alpha_imag=imag(alpha),
        overlap=overlap,backward_error=first(result.backward_errors),
        target_real=real(target),target_imag=imag(target),
        evidence_scope=String(case.evidence_scope),
    )
    @printf("%-8s %-12s %2d N=%d R=%8.4f alpha=%.10f%+.10fi overlap=%s berr=%.2e\n",
        flow.flow,stage,index,N,R,real(alpha),imag(alpha),
        isnan(overlap) ? "seed" : @sprintf("%.6f",overlap),
        first(result.backward_errors))
    return (row=row,case=case,alpha=alpha,
            vector=copy(view(result.vectors,:,1)))
end

function mach_seed(flow,rows;N=49)
    target = ComplexF64(flow.alpha_ref)
    vector = nothing
    point = nothing
    for (index,Mr) in enumerate(MACH_PATH)
        point = solve_point(flow;Mr=Mr,R=flow.R,N=N,target=target,
            previous_vector=vector,stage="mach",index=index)
        push!(rows,point.row)
        write_rows(joinpath(OUTDIR,"continuation_points.csv"),rows)
        target,vector = point.alpha,point.vector
    end
    return point
end

function radial_path(flow,start,R_values,stage,rows)
    points = [start]
    target,vector = start.alpha,start.vector
    for (index,R) in enumerate(R_values)
        point = solve_point(flow;Mr=0.30,R=R,N=49,target=target,
            previous_vector=vector,stage=stage,index=index)
        push!(points,point)
        push!(rows,point.row)
        write_rows(joinpath(OUTDIR,"continuation_points.csv"),rows)
        target,vector = point.alpha,point.vector
    end
    return points
end

function sign_brackets(points)
    brackets = Tuple{Any,Any}[]
    for index in 1:length(points)-1
        a,b = points[index],points[index+1]
        signbit(imag(a.alpha)) != signbit(imag(b.alpha)) && push!(brackets,(a,b))
    end
    return brackets
end

function refine_bracket(flow,left,right,rows,neutral_rows;max_iterations=18)
    imag(left.alpha)*imag(right.alpha) <= 0 || error("invalid neutral bracket")
    left.row.R > right.row.R && ((left,right)=(right,left))
    point = abs(imag(left.alpha)) < abs(imag(right.alpha)) ? left : right
    for index in 1:max_iterations
        fleft,fright = imag(left.alpha),imag(right.alpha)
        trial_R = (left.row.R*fright-right.row.R*fleft)/(fright-fleft)
        width = right.row.R-left.row.R
        if !isfinite(trial_R) || trial_R <= left.row.R+0.1width ||
           trial_R >= right.row.R-0.1width
            trial_R = (left.row.R+right.row.R)/2
        end
        seed = abs(trial_R-left.row.R) <= abs(trial_R-right.row.R) ? left : right
        point = solve_point(flow;Mr=0.30,R=trial_R,N=49,target=seed.alpha,
            previous_vector=seed.vector,stage="refine",index=index)
        push!(rows,point.row)
        write_rows(joinpath(OUTDIR,"continuation_points.csv"),rows)
        if abs(imag(point.alpha)) < 1e-8 || width < 1e-7
            break
        elseif signbit(imag(left.alpha)) != signbit(imag(point.alpha))
            right = point
        else
            left = point
        end
    end

    for N in (49,59,69)
        check = N == 49 ? point : solve_point(flow;Mr=0.30,R=point.row.R,N=N,
            target=point.alpha,previous_vector=nothing,stage="grid_check",index=N)
        push!(neutral_rows,(
            flow=flow.flow,Ro=flow.Ro,Mr=0.30,N=N,R=point.row.R,
            beta=flow.beta,alpha_real=real(check.alpha),alpha_imag=imag(check.alpha),
            backward_error=check.row.backward_error,
            fixed_beta_neutral_at_N49=abs(imag(point.alpha)) < 1e-7,
            evidence_scope=check.row.evidence_scope,
        ))
        write_rows(joinpath(OUTDIR,"neutral_candidates.csv"),neutral_rows)
    end
end

mkpath(OUTDIR)
rows = NamedTuple[]
neutral_rows = NamedTuple[]
summary_rows = NamedTuple[]
for flow in FLOWS
    seed = mach_seed(flow,rows)
    lower = radial_path(flow,seed,flow.lower,"R_down",rows)
    upper = radial_path(flow,seed,flow.upper,"R_up",rows)
    brackets = vcat(sign_brackets(lower),sign_brackets(upper))
    push!(summary_rows,(
        flow=flow.flow,Ro=flow.Ro,Mr=0.30,N=49,beta=flow.beta,
        R_min=minimum(vcat(flow.lower,[flow.R],flow.upper)),
        R_max=maximum(vcat(flow.lower,[flow.R],flow.upper)),
        sign_change_count=length(brackets),
        reference_alpha_real=real(seed.alpha),reference_alpha_imag=imag(seed.alpha),
    ))
    write_rows(joinpath(OUTDIR,"scan_summary.csv"),summary_rows)
    for (left,right) in brackets
        refine_bracket(flow,left,right,rows,neutral_rows)
    end
end

open(joinpath(OUTDIR,"metadata.txt"),"w") do io
    println(io,"study=corrected_compressible_BEK_fixed_beta_neutral_continuation")
    println(io,"date=2026-09-12")
    println(io,"neutral_condition=imag(alpha)=0 at fixed beta")
    println(io,"Mr_definition=local_differential_rotation_Mach_r_absDeltaOmega_over_ainf")
    println(io,"operator_source=$(FRONTEND)")
    println(io,"legacy_backend_preserved=$(BACKEND)")
end
println("Wrote $(OUTDIR)")

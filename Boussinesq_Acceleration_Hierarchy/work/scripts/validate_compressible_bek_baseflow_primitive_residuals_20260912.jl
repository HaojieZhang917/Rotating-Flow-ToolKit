using BSplineKit
using LinearAlgebra
using Printf
using PyCall

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const TOOLKIT_ROOT=normpath(joinpath(ROOT,".."))
const BACKEND=joinpath(TOOLKIT_ROOT,"RotatingDiskFlow","src","CRD_STA.jl")
const FRONTEND=joinpath(ROOT,"work","src","CompressibleBEKLinearStability.jl")
const RESIDUAL_SOURCE=joinpath(ROOT,"work","src","CompressibleBEKBaseflowResiduals.jl")
const OUTDIR=joinpath(ROOT,"work","results",
    "compressible_bek_baseflow_primitive_residual_gate_v3_20260915")

function load_solver()
    backend_text=replace(read(BACKEND,String),
        "using DifferentialEquations"=>"",
        "using BoundaryValueDiffEq"=>"")
    include_string(Main,backend_text,BACKEND*"#baseflow-residual-gate")
    frontend_text=replace(read(FRONTEND,String),
        "using DifferentialEquations"=>"",
        "using BoundaryValueDiffEq"=>"")
    old_block="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at \$BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new_block="""const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const CRD = Main"""
    occursin(old_block,frontend_text) || error("frontend loader block changed")
    include_string(Main,replace(frontend_text,old_block=>new_block),
        FRONTEND*"#baseflow-residual-gate")
    Main.CompressibleBEKLinearStability
end

const CBEK=load_solver()
include(RESIDUAL_SOURCE)
const BRES=CompressibleBEKBaseflowResiduals
const FLOWS=[
    (flow="von_Karman",Ro=-1.0,R=285.36),
    (flow="Ekman",Ro=0.0,R=116.25),
    (flow="Bodewadt",Ro=1.0,R=27.26),
]
const THERMAL_CASES=[
    (label="viscous_heating",Mr=0.3,Tw=1.0),
    (label="wall_heat_transfer",Mr=0.0,Tw=0.8),
]
const GRIDS=(80,120,160,200)
const REGIONS=((label="core",max_y=12.0),(label="extended",max_y=39.0))

function sampled_thermal_functions(source,y)
    f_spline=BSplineKit.interpolate(source.y,source.f,BSplineOrder(4))
    q_spline=BSplineKit.interpolate(source.y,source.q,BSplineOrder(4))
    f=[point<=last(source.y) ? f_spline(point) : 0.0 for point in y]
    q=[point<=last(source.y) ? q_spline(point) : 0.0 for point in y]
    f[1]=0.0; f[end]=0.0; q[1]=1.0; q[end]=0.0
    return f,q
end

function write_rows(path,rows)
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(row->println(io,join(values(row),',')),rows)
    end
end

mkpath(OUTDIR)
rows=NamedTuple[]
profile_rows=NamedTuple[]
residual_profile_rows=NamedTuple[]
for flow in FLOWS, thermal in THERMAL_CASES, N in GRIDS
    p=CBEK.PhysicalParameters(Mr=thermal.Mr,Tw=thermal.Tw,Ro=flow.Ro,
        R=flow.R,beta=0.0,omega=0.0,gamma=1.4,Pr=0.72)
    n=CBEK.NumericalParameters(N=N,domain=:finite,ymax=40.0,neigs=1)
    case=CBEK.prepare_case(p,n;allow_unverified_ro=true)
    source=CBEK.source_profiles(p.Ro,p.Pr)
    f,q=sampled_thermal_functions(source,case.y)
    checks=NamedTuple[]
    for region in REGIONS
        check=BRES.baseflow_primitive_residuals(
            case.y,case.D,case.D2,case.F,case.G,case.H./case.T,
            case.T,case.rho,f,q;Ro=p.Ro,Mr=p.Mr,Tw=p.Tw,
            gamma=p.gamma,Pr=p.Pr,max_y=region.max_y)
        push!(checks,check)
        for name in keys(check.metrics)
            metric=getfield(check.metrics,name)
            push!(rows,(
                flow=flow.flow,thermal_case=thermal.label,Ro=flow.Ro,Co=check.Co,
                Mr=thermal.Mr,Tw=thermal.Tw,N=N,region=region.label,
                equation=String(name),absolute_linf=metric.absolute_linf,
                relative_linf=metric.relative_linf,rms=metric.rms,
                term_relative_linf=metric.term_relative_linf,
                scale_linf=metric.scale_linf,max_y=check.max_y,
            ))
        end
    end
    q_spline=BSplineKit.interpolate(case.y,q,BSplineOrder(4))
    q_values=map(point->q_spline(point),(1.0,5.0,10.0,20.0,30.0,35.0,39.0))
    half_index=argmin(abs.(q .- 0.5))
    push!(profile_rows,(
        flow=flow.flow,thermal_case=thermal.label,Ro=flow.Ro,Mr=thermal.Mr,
        Tw=thermal.Tw,N=N,min_T=minimum(case.T),max_T=maximum(case.T),
        max_state_error=maximum(abs.(case.rho.*case.T.-1)),
        thermal_bc_error=case.thermal_bc_error,q_y1=q_values[1],q_y5=q_values[2],
        q_y10=q_values[3],q_y20=q_values[4],q_y30=q_values[5],
        q_y35=q_values[6],q_y39=q_values[7],q_half_y=case.y[half_index],
    ))
    core,extended=checks
    if N==120
        for index in eachindex(case.y)
            r=extended.residuals
            push!(residual_profile_rows,(
                flow=flow.flow,thermal_case=thermal.label,Ro=flow.Ro,
                Mr=thermal.Mr,Tw=thermal.Tw,N=N,index=index,y=case.y[index],
                U=case.F[index],V=case.G[index],W=(case.H./case.T)[index],
                T=case.T[index],rho=case.rho[index],
                continuity=r.primitive_continuity[index],
                radial_momentum=r.primitive_radial[index],
                azimuthal_momentum=r.primitive_azimuthal[index],
                state=r.primitive_state[index],energy=r.primitive_energy[index],
                current_energy=r.current_energy[index],
                literature_f=r.literature_f[index],literature_q=r.literature_q[index],
                current_f=r.current_f[index],current_q=r.current_q[index],
            ))
        end
    end
    @printf("%-11s %-18s N=%3d vel=(%.2e,%.2e,%.2e) thermal(lit,current)=(%.2e,%.2e) extended-lit=%.2e\n",
        flow.flow,thermal.label,N,
        core.metrics.primitive_continuity.relative_linf,
        core.metrics.primitive_radial.relative_linf,
        core.metrics.primitive_azimuthal.relative_linf,
        core.metrics.primitive_energy.relative_linf,
        core.metrics.current_energy.relative_linf,
        extended.metrics.primitive_energy.relative_linf)
    write_rows(joinpath(OUTDIR,"residual_metrics.csv"),rows)
    write_rows(joinpath(OUTDIR,"profile_checks.csv"),profile_rows)
    !isempty(residual_profile_rows) && write_rows(
        joinpath(OUTDIR,"residual_profiles_N120.csv"),residual_profile_rows)
end
open(joinpath(OUTDIR,"metadata.txt"),"w") do io
    println(io,"study=compressible_BEK_baseflow_pre_DH_primitive_residual_gate")
    println(io,"date=2026-09-12")
    println(io,"coordinate_for_residual=physical_normal_eta_via_deta_equals_rho_Dy")
    println(io,"evaluation_regions=interior_nodes_y_le_12_and_y_le_39")
    println(io,"relative_norm=absolute_linf/max(term_scale_linf,1)")
    println(io,"term_relative_norm=absolute_linf/max(term_scale_linf,eps(Float64))")
    println(io,"pointwise_residual_profiles=N120_all_nodes")
    println(io,"operator_source=$(FRONTEND)")
    println(io,"residual_source=$(RESIDUAL_SOURCE)")
    println(io,"legacy_backend_preserved=$(BACKEND)")
end
println("Wrote $(OUTDIR)")

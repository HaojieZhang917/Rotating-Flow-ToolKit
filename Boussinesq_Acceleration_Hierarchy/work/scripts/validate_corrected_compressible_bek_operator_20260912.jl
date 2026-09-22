using BoundaryValueDiffEq
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
const OUTDIR = joinpath(ROOT, "work", "results", "compressible_bek_corrected_operator_v8_20260912")

function load_solver()
    backend_text = replace(read(BACKEND, String),
        "using DifferentialEquations" => "using BoundaryValueDiffEq")
    include_string(Main, backend_text, BACKEND * "#minimal-import-corrected-operator")

    frontend_text = replace(read(FRONTEND, String),
        "using DifferentialEquations" => "using BoundaryValueDiffEq")
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
        FRONTEND * "#minimal-import-corrected-operator")
    return Main.CompressibleBEKLinearStability
end

const CBEK = load_solver()
const REFERENCES = [
    (flow="von_Karman", Ro=-1.0, R=285.36, beta=0.07759, alpha=0.38482),
    (flow="Ekman", Ro=0.0, R=116.25, beta=0.13200, alpha=0.27750),
    (flow="Bodewadt", Ro=1.0, R=27.26, beta=0.12700, alpha=0.52310),
]

function make_case(ref; Mr=0.01, R=ref.R, target=ComplexF64(ref.alpha), N=49,
                   neigs=1)
    p = CBEK.PhysicalParameters(
        Mr=Mr, Tw=1.0, Ro=ref.Ro, R=R, beta=ref.beta,
        omega=0.0, gamma=1.4, Pr=0.72,
    )
    n = CBEK.NumericalParameters(
        N=N, domain=:rational, target=target, neigs=neigs,
        tol=1e-9, maxit=500, balance=true, balance_iterations=2,
        regularization=0.0,
    )
    CBEK.prepare_case(p, n; allow_unverified_ro=true)
end

function fixed_reference_audit(; mach_values=(0.01,0.30), grids=(39,49,59))
    rows = NamedTuple[]
    for ref in REFERENCES, Mr in mach_values, N in grids
        case = make_case(ref; Mr=Mr, N=N)
        gates = CBEK.operator_structure_gates(case)
        result = CBEK.solve_mode(case)
        alpha = first(result.alpha)
        push!(rows, (
            flow=ref.flow, Ro=ref.Ro, Mr=Mr, N=N, R=ref.R,
            beta=ref.beta, alpha_ref=ref.alpha,
            alpha_real=real(alpha), alpha_imag=imag(alpha),
            backward_error=first(result.backward_errors),
            continuity_gate=gates.continuity_curvature,
            continuity_transport_gate=gates.continuity_base_transport,
            radial_coriolis_gate=gates.radial_coriolis,
            azimuthal_coriolis_gate=gates.azimuthal_coriolis,
            energy_gate=gates.energy_wallnormal_scale,
            direct_Co_gate=gates.bodewadt_direct_Co,
        ))
        @printf("fixed %-11s Mr=%.2f N=%d alpha=%.11f%+.11fi berr=%.2e\n",
            ref.flow,Mr,N,real(alpha),imag(alpha),first(result.backward_errors))
    end
    return rows
end

function write_rows(path,rows)
    open(path,"w") do io
        println(io,join(keys(first(rows)),','))
        foreach(row -> println(io,join(values(row),',')),rows)
    end
end

mkpath(OUTDIR)
rows = fixed_reference_audit()
write_rows(joinpath(OUTDIR,"fixed_reference_roots.csv"),rows)
open(joinpath(OUTDIR,"metadata.txt"),"w") do io
    println(io,"study=corrected_compressible_BEK_operator_fixed_reference_audit")
    println(io,"date=2026-09-12")
    println(io,"Mr_definition=local_differential_rotation_Mach_r_absDeltaOmega_over_ainf")
    println(io,"operator_source=$(FRONTEND)")
    println(io,"legacy_backend_preserved=$(BACKEND)")
end
println("Wrote $(OUTDIR)")

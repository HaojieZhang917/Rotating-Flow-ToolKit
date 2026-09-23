using BoundaryValueDiffEq
using BSplineKit
using LinearAlgebra
using NonlinearEigenproblems
using Printf
using PyCall

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TOOLKIT_ROOT = normpath(joinpath(ROOT, ".."))
const BACKEND = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "CRD_STA.jl")
const FRONTEND = joinpath(ROOT, "work", "src", "CompressibleBEKLinearStability.jl")
const OUTDIR = joinpath(ROOT, "work", "results", "ro0_matched_notebook_comparison_20260922")

function load_solver()
    backend_text = replace(read(BACKEND, String),
        "using DifferentialEquations" => "using BoundaryValueDiffEq")
    include_string(Main, backend_text, BACKEND * "#matched-ro0")
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
    include_string(Main, replace(frontend_text, old_block => new_block),
        FRONTEND * "#matched-ro0")
    Main.CompressibleBEKLinearStability
end

const CBEK = load_solver()

function main()
    mkpath(OUTDIR)
    p = CBEK.PhysicalParameters(
        Ro=0.0, Tw=1.0, Mr=0.1, R=116.3, beta=0.137,
        omega=0.0, gamma=1.4, Pr=0.72,
    )
    n = CBEK.NumericalParameters(
        N=69, domain=:rational, ymax=40.0, target=0.54+0im,
        neigs=6, tol=1e-11, maxit=800, balance=true,
        balance_iterations=2, regularization=0.0,
    )
    c = CBEK.prepare_case(p, n; allow_unverified_ro=true)
    result = CBEK.solve_mode(c)
    open(joinpath(OUTDIR, "modes.csv"), "w") do io
        println(io,"mode,alpha_real,alpha_imag,backward_error,distance_to_legacy_root")
        for j in eachindex(result.alpha)
            alpha = result.alpha[j]
            @printf("%d %.12f%+.12fi residual=%.3e\n",j,real(alpha),imag(alpha),result.backward_errors[j])
            println(io,join((j,real(alpha),imag(alpha),result.backward_errors[j],
                abs(alpha-(0.5403612112632648+0.0003372888856162062im))),','))
        end
    end
    @printf("max|T-1|=%.8e evidence_scope=%s\n",
        maximum(abs.(c.T .- 1)),string(c.evidence_scope))
end

main()

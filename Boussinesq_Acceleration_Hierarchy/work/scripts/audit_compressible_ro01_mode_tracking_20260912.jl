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
const OUTDIR = joinpath(ROOT, "work", "results", "compressible_solver_ro01_audit_20260912")

"""Load the existing operator without the broken DifferentialEquations meta-package.

The replacement changes imports only.  The base-flow, interpolation, operator,
boundary elimination, balancing and eigensolver function bodies are otherwise
evaluated verbatim from the registered source files.
"""
function load_registered_solver()
    backend_text = read(BACKEND, String)
    backend_text = replace(backend_text,
        "using DifferentialEquations" => "using BoundaryValueDiffEq")
    include_string(Main, backend_text, BACKEND * "#minimal-import-audit")

    frontend_text = read(FRONTEND, String)
    frontend_text = replace(frontend_text,
        "using DifferentialEquations" => "using BoundaryValueDiffEq")
    old_block = """const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\", \"..\"))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, \"RotatingDiskFlow\", \"src\", \"RotatingDiskFlow.jl\")
isfile(BACKEND_SOURCE) || error(\"RotatingDiskFlow backend not found at \$BACKEND_SOURCE\")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner"""
    new_block = """const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, \"..\", \"..\", \"..\"))
const CRD = Main"""
    occursin(old_block, frontend_text) || error("frontend backend block changed; audit loader must be reviewed")
    frontend_text = replace(frontend_text, old_block => new_block)
    include_string(Main, frontend_text, FRONTEND * "#minimal-import-audit")
    return Main.CompressibleBEKLinearStability
end

const CBEK = load_registered_solver()

const FLOWS = [
    (flow="Ekman", Ro=0.0, R=116.25, beta=0.132, alpha_ref=0.2775),
    (flow="Bodewadt", Ro=1.0, R=27.26, beta=0.127, alpha_ref=0.5231),
]
const MACH_PATH = [0.01, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
const N = 69

function grid_mapping_audit(; grids=(59, 69, 79, 89, 99))
    rows = NamedTuple[]
    for grid_N in grids
        numerical = CBEK.NumericalParameters(N=grid_N, domain=:rational)
        D, D2, y = CBEK.stability_grid(numerical)
        legacy_D, legacy_D2, legacy_y = Main.CRD_BF.Cheb(
            grid_N; domain=:rational, ymax=40.0,
        )
        finite_y = y[isfinite.(y)]
        push!(rows, (
            N=grid_N,
            node_count=length(y),
            finite_node_count=length(finite_y),
            finite_nodes_strictly_increasing=all(diff(finite_y) .> 0),
            infinity_endpoint=isinf(last(y)),
            last_finite_y=last(finite_y),
            finite_nodes_beyond_BVP_support=count(>(40.0), finite_y),
            legacy_repeated_y40=count(==(40.0), vec(legacy_y)),
            relative_D_change=norm(D-legacy_D) / norm(legacy_D),
            relative_D2_change=norm(D2-legacy_D2) / norm(legacy_D2),
        ))
    end
    open(joinpath(OUTDIR, "semi_infinite_grid_mapping.csv"), "w") do io
        println(io, join(keys(first(rows)), ','))
        foreach(row -> println(io, join(values(row), ',')), rows)
    end
    return rows
end

function parameters(flow, Mr; target, neigs=1, grid_N=N, domain=:rational)
    physical = CBEK.PhysicalParameters(
        Mr=Mr, Tw=1.0, Ro=flow.Ro, R=flow.R,
        beta=flow.beta, omega=0.0, gamma=1.4, Pr=0.72,
    )
    numerical = CBEK.NumericalParameters(
        N=grid_N, domain=domain, ymax=40.0,
        target=ComplexF64(target), neigs=neigs, tol=1e-10,
        maxit=800, balance=true, balance_iterations=2,
        regularization=0.0,
    )
    return physical, numerical
end

function endpoint_grid_audit(; grids=(59, 69, 79), flows=FLOWS,
                             filename="endpoint_grid_convergence_Mr03.csv")
    tracking_path = joinpath(OUTDIR, "mach_mode_tracking_N69.csv")
    isfile(tracking_path) || error("run the Mach tracking study before the endpoint grid audit")
    data, header = readdlm(tracking_path, ',', Any, '\n'; header=true)
    names = vec(header)
    flow_column = findfirst(==("flow"), names)
    mach_column = findfirst(==("Mr"), names)
    real_column = findfirst(==("alpha_real"), names)
    imag_column = findfirst(==("alpha_imag"), names)
    rows = NamedTuple[]
    for flow in flows
        selected = findfirst(eachindex(axes(data, 1))) do row
            data[row, flow_column] == flow.flow && data[row, mach_column] == 0.3
        end
        isnothing(selected) && error("missing Mr=0.3 tracking endpoint for $(flow.flow)")
        target = ComplexF64(data[selected, real_column], data[selected, imag_column])
        for grid_N in grids
            p, n = parameters(flow, 0.3; target=target, grid_N=grid_N)
            case = CBEK.prepare_case(p, n; allow_unverified_ro=true)
            result = CBEK.solve_mode(case)
            alpha = first(result.alpha)
            push!(rows, (
                flow=flow.flow, Ro=flow.Ro, N=grid_N, domain="rational",
                alpha_real=real(alpha), alpha_imag=imag(alpha),
                backward_error=first(result.backward_errors),
                distance_to_N69_tracked_target=abs(alpha-target),
                evidence_scope=String(case.evidence_scope),
            ))
            @printf("grid %s N=%d alpha=%.12f%+.12fi berr=%.3e\n",
                flow.flow, grid_N, real(alpha), imag(alpha), first(result.backward_errors))
        end
    end
    open(joinpath(OUTDIR, filename), "w") do io
        println(io, join(keys(first(rows)), ','))
        foreach(row -> println(io, join(values(row), ',')), rows)
    end
    return rows
end

function retained_indices(degree)
    n = degree + 1
    removed = (n, n+1, 2n, 2n+1, 3n, 3n+1, 4n, 4n+1, 5n)
    return setdiff(1:5n, removed)
end

function interpolate_reduced_mode(old_case, old_vector, new_case)
    old_n = old_case.numerical.N + 1
    new_n = new_case.numerical.N + 1
    old_full = zeros(ComplexF64, 5old_n)
    old_full[retained_indices(old_case.numerical.N)] .= old_vector
    new_full = zeros(ComplexF64, 5new_n)
    # Interpolate in the compact coordinate s=y/(1+y), with infinity mapped
    # exactly to s=1.  This avoids extrapolation between the last finite node
    # and the asymptotic boundary when changing N.
    compact(point) = isinf(point) ? 1.0 : point / (1 + point)
    sample_y = compact.(old_case.y)
    new_y = compact.(new_case.y)
    function linear_sample(values, point)
        point <= first(sample_y) && return first(values)
        point >= last(sample_y) && return last(values)
        left = searchsortedlast(sample_y, point)
        fraction = (point-sample_y[left]) / (sample_y[left+1]-sample_y[left])
        return (1-fraction)*values[left] + fraction*values[left+1]
    end
    for field in 0:4
        old_block = view(old_full, field*old_n+1:(field+1)*old_n)
        new_full[field*new_n+1:(field+1)*new_n] .=
            [linear_sample(old_block, point) for point in new_y]
    end
    reduced = new_full[retained_indices(new_case.numerical.N)]
    reduced ./= norm(reduced)
    return reduced
end

function bodewadt_grid_overlap_continuation()
    flow = last(FLOWS)
    grids = (59, 69, 79, 89, 99)
    target = ComplexF64(0.5634677931234605, -0.15051884586747033)
    previous_case = nothing
    previous_vector = nothing
    rows = NamedTuple[]
    for grid_N in grids
        p, n = parameters(flow, 0.3; target=target, grid_N=grid_N)
        case = CBEK.prepare_case(p, n; allow_unverified_ro=true)
        initial = isnothing(previous_case) ? nothing :
            interpolate_reduced_mode(previous_case, previous_vector, case)
        result = CBEK.solve_mode(case; initial_vector=initial)
        alpha = first(result.alpha)
        overlap = first(result.overlaps)
        push!(rows, (
            flow=flow.flow, Ro=flow.Ro, N=grid_N, domain="rational",
            alpha_real=real(alpha), alpha_imag=imag(alpha), overlap=overlap,
            backward_error=first(result.backward_errors),
            step_change=abs(alpha-target), evidence_scope=String(case.evidence_scope),
        ))
        @printf("grid-track %s N=%d alpha=%.12f%+.12fi overlap=%s berr=%.3e\n",
            flow.flow, grid_N, real(alpha), imag(alpha),
            isnan(overlap) ? "seed" : @sprintf("%.6f", overlap),
            first(result.backward_errors))
        target = alpha
        previous_case = case
        previous_vector = copy(view(result.vectors, :, 1))
    end
    open(joinpath(OUTDIR, "bodewadt_grid_overlap_continuation_Mr03.csv"), "w") do io
        println(io, join(keys(first(rows)), ','))
        foreach(row -> println(io, join(values(row), ',')), rows)
    end
    return rows
end

function thermal_audit()
    rows = NamedTuple[]
    for flow in FLOWS
        p, n = parameters(flow, 0.3; target=flow.alpha_ref, neigs=1)
        case = CBEK.prepare_case(p, n; allow_unverified_ro=true)
        checks = CBEK.validate_result(CBEK.StabilityResult(
            case, ComplexF64[], zeros(ComplexF64, 5N - 4, 0), Float64[], Float64[],
        ))
        push!(rows, (
            flow=flow.flow, Ro=flow.Ro, min_T=minimum(case.T),
            max_T=maximum(case.T), max_abs_T_minus_one=maximum(abs.(case.T .- 1)),
            max_abs_rho_minus_one=maximum(abs.(case.rho .- 1)),
            thermal_bc_error=case.thermal_bc_error,
            wall_T_error=checks.base_bc.wall_T, far_T_error=checks.base_bc.far_T,
        ))
    end
    open(joinpath(OUTDIR, "corrected_thermal_baseflow.csv"), "w") do io
        println(io, join(keys(first(rows)), ','))
        foreach(row -> println(io, join(values(row), ',')), rows)
    end
    return rows
end

function track_flow(flow)
    rows = NamedTuple[]
    previous_vector = nothing
    target = ComplexF64(flow.alpha_ref)
    for (step_index, Mr) in enumerate(MACH_PATH)
        p, n = parameters(flow, Mr; target=target)
        case = CBEK.prepare_case(p, n; allow_unverified_ro=true)
        result = CBEK.solve_mode(case; initial_vector=previous_vector)
        alpha = first(result.alpha)
        overlap = first(result.overlaps)
        backward_error = first(result.backward_errors)
        push!(rows, (
            flow=flow.flow, Ro=flow.Ro, step=step_index, Mr=Mr,
            Ma=case.Ma, target_real=real(target), target_imag=imag(target),
            alpha_real=real(alpha), alpha_imag=imag(alpha),
            overlap=overlap, backward_error=backward_error,
            distance_to_supplied_reference=abs(alpha-flow.alpha_ref),
            neutral_residual=imag(alpha), evidence_scope=String(case.evidence_scope),
        ))
        @printf("%s Mr=%.2f alpha=%.12f%+.12fi overlap=%s berr=%.3e\n",
            flow.flow, Mr, real(alpha), imag(alpha),
            isnan(overlap) ? "seed" : @sprintf("%.6f", overlap), backward_error)
        previous_vector = copy(view(result.vectors, :, 1))
        target = alpha
    end
    return rows
end

function direct_ro_branch_audit()
    p = CBEK.PhysicalParameters(
        Mr=0.1, Tw=1.0, Ro=-1.0, R=285.36,
        beta=0.07759, omega=0.0, gamma=1.4, Pr=0.72,
    )
    n = CBEK.NumericalParameters(N=49, target=0.38482+0im, neigs=1,
        regularization=0.0)
    case = CBEK.prepare_case(p, n; allow_unverified_ro=true)
    unified = CBEK.assemble_qep(case)
    direct_at = ro -> begin
        coefficients = Main.Spatial_mode_BEK(
            case.F, case.G, case.H, case.rho, case.lambda, case.kappa, case.T,
            p.Pr, p.gamma, p.R, case.Ma, n.N, ro, case.Co, case.D, case.D2;
            regularization=0.0,
        )
        raw = Main.assemble_mat(coefficients, case.D, case.D2, p.beta, p.omega)
        Main.boundary_condition(raw..., n.N)
    end
    direct_minus_one = direct_at(-1.0)
    direct_left_limit = direct_at(-1.0 + 1e-8)
    relative(a, b) = norm(a-b) / max(norm(b), eps(Float64))
    rows = [(matrix="A$(j-1)",
             direct_special_vs_left_limit=relative(direct_minus_one[j], direct_left_limit[j]),
             unified_vs_direct_special=relative(unified[j], direct_minus_one[j]),
             unified_vs_left_limit=relative(unified[j], direct_left_limit[j])) for j in 1:3]
    open(joinpath(OUTDIR, "ro_special_branch_discontinuity.csv"), "w") do io
        println(io, join(keys(first(rows)), ','))
        foreach(row -> println(io, join(values(row), ',')), rows)
    end
    return rows
end

mkpath(OUTDIR)
if "--grid-map-only" in ARGS
    grid_mapping_audit()
elseif "--grid-only" in ARGS
    endpoint_grid_audit()
elseif "--bodewadt-high-grid-only" in ARGS
    endpoint_grid_audit(; grids=(89, 99), flows=(last(FLOWS),),
        filename="endpoint_grid_convergence_Mr03_Bodewadt_N89_N99.csv")
elseif "--bodewadt-grid-track-only" in ARGS
    bodewadt_grid_overlap_continuation()
else
    open(joinpath(OUTDIR, "runtime_metadata.txt"), "w") do io
        println(io, "study=compressible_Ro01_mode_tracking_audit")
        println(io, "date=2026-09-12")
        println(io, "julia_version=$(VERSION)")
        println(io, "backend_source=$(BACKEND)")
        println(io, "frontend_source=$(FRONTEND)")
        println(io, "runtime_workaround=replace_using_DifferentialEquations_with_using_BoundaryValueDiffEq_in_memory")
        println(io, "original_sources_modified=false")
        println(io, "N=$(N)")
        println(io, "mach_path=$(join(MACH_PATH, ';'))")
    end

    thermal_audit()
    direct_ro_branch_audit()
    tracking_rows = reduce(vcat, track_flow(flow) for flow in FLOWS)
    open(joinpath(OUTDIR, "mach_mode_tracking_N69.csv"), "w") do io
        println(io, join(keys(first(tracking_rows)), ','))
        foreach(row -> println(io, join(values(row), ',')), tracking_rows)
    end
end

println("Wrote $(OUTDIR)")

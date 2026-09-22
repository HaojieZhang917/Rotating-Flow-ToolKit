using DelimitedFiles
using Printf

include(joinpath(@__DIR__, "..", "src", "CompressibleBEKLinearStability.jl"))
using .CompressibleBEKLinearStability

const OUTDIR = normpath(joinpath(
    @__DIR__, "..", "results", "compressible_bek_lsa_validation_20260911",
))

p = PhysicalParameters(
    Mr=0.1, Tw=1.0, Ro=-1.0, R=285.36,
    beta=0.07759, omega=0.0, gamma=1.4, Pr=0.72,
)
grids = [
    (domain=:rational, N=59, ymax=40.0),
    (domain=:rational, N=69, ymax=40.0),
    (domain=:rational, N=79, ymax=40.0),
    (domain=:finite, N=79, ymax=40.0),
    (domain=:finite, N=99, ymax=40.0),
]

rows = NamedTuple[]
for grid in grids
    n = NumericalParameters(
        N=grid.N, domain=grid.domain, ymax=grid.ymax,
        target=0.38482 + 0im, neigs=1, tol=1e-11,
        balance=true, regularization=0.0,
    )
    @printf("Solving domain=%s N=%d\n", String(grid.domain), grid.N)
    c = prepare_case(p, n)
    audit = audit_expanded_operator(c)
    result = solve_mode(c)
    checks = validate_result(result)
    row = (
        domain=grid.domain, N=grid.N, ymax=grid.ymax,
        alpha_real=real(result.alpha[1]), alpha_imag=imag(result.alpha[1]),
        backward_error=result.backward_errors[1],
        operator_A0=audit.relative[1], operator_A1=audit.relative[2],
        operator_A2=audit.relative[3],
        thermal_bc_error=checks.thermal_bc_error,
        wall_bc_max=maximum(values(checks.base_bc)[1:4]),
        far_bc_max=maximum(values(checks.base_bc)[5:7]),
        Tmin=minimum(c.T), Tmax=maximum(c.T),
    )
    push!(rows, row)
    @printf("  alpha=%.15g%+.15gi, berr=%.3e, operator max=%.3e\n",
            row.alpha_real, row.alpha_imag, row.backward_error,
            maximum((row.operator_A0, row.operator_A1, row.operator_A2)))
end

mkpath(OUTDIR)
open(joinpath(OUTDIR, "grid_convergence.csv"), "w") do io
    println(io, join(keys(rows[1]), ','))
    for row in rows
        println(io, join(values(row), ','))
    end
end
open(joinpath(OUTDIR, "README.txt"), "w") do io
    println(io, "Strict-paper von Karman near-neutral benchmark")
    println(io, "Mr=0.1, Tw=1.0, Ro=-1, R=285.36, beta=0.07759, omega=0")
    println(io, "Ma_internal=Mr/R=$(p.Mr/p.R)")
    println(io, "published rounded target alpha_r=0.38482; exact neutrality is not assumed")
    println(io, "all profiles sampled on Dorodnitsyn-Howarth similarity coordinate y")
    println(io, "all QEPs use zero artificial regularization")
end

println("Wrote $OUTDIR")

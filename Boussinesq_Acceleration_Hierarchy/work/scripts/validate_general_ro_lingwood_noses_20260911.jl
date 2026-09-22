using DelimitedFiles
using Printf

include(joinpath(@__DIR__, "..", "src", "CompressibleBEKLinearStability.jl"))
using .CompressibleBEKLinearStability

const OUTDIR = normpath(joinpath(
    @__DIR__, "..", "results", "compressible_bek_general_ro_lingwood_20260911",
))

# Approximate stationary-wave marginal-curve noses from Lingwood (1997),
# table 2. These are graphical/rounded reference values, not exact roots.
references = [
    (Ro=-1.0, R=290.1, alpha=0.381, beta=0.077),
    (Ro=-0.5, R=160.9, alpha=0.467, beta=0.117),
    (Ro= 0.0, R=116.3, alpha=0.528, beta=0.137),
    (Ro= 0.5, R= 75.9, alpha=0.544, beta=0.140),
    (Ro= 1.0, R= 27.4, alpha=0.487, beta=0.115),
]

rows = NamedTuple[]
for ref in references
    p = PhysicalParameters(
        Mr=0.1, Tw=1.0, Ro=ref.Ro, R=ref.R,
        beta=ref.beta, omega=0.0, gamma=1.4, Pr=0.72,
    )
    n = NumericalParameters(
        N=69, domain=:rational, ymax=40.0,
        target=ComplexF64(ref.alpha), neigs=2, tol=1e-11,
        regularization=0.0, balance=true,
    )
    @printf("Ro=%+.1f R=%.1f beta=%.3f target=%.3f\n",
            ref.Ro, ref.R, ref.beta, ref.alpha)
    case = prepare_case(p, n; allow_unverified_ro=true)
    audit = audit_expanded_operator(case)
    result = solve_mode(case)
    j = argmin(abs.(result.alpha .- ref.alpha))
    row = (
        Ro=ref.Ro, R_reference=ref.R, beta_reference=ref.beta,
        alpha_reference=ref.alpha,
        alpha_real=real(result.alpha[j]), alpha_imag=imag(result.alpha[j]),
        alpha_relative_difference=abs(real(result.alpha[j])-ref.alpha)/ref.alpha,
        backward_error=result.backward_errors[j],
        operator_relative_max=maximum(audit.relative),
        evidence_scope=case.evidence_scope,
    )
    push!(rows, row)
    @printf("  alpha=%.12g%+.12gi relative(real)=%.3e backward=%.3e\n",
            row.alpha_real, row.alpha_imag, row.alpha_relative_difference,
            row.backward_error)
end

mkpath(OUTDIR)
open(joinpath(OUTDIR, "lingwood_nose_diagnostic.csv"), "w") do io
    println(io, join(keys(rows[1]), ','))
    foreach(row -> println(io, join(values(row), ',')), rows)
end
open(joinpath(OUTDIR, "README.txt"), "w") do io
    println(io, "Low-Mach diagnostic against approximate Lingwood (1997) table-2 noses")
    println(io, "Mr=0.1, Tw=1, omega=0, Pr=0.72, gamma=1.4, N=69")
    println(io, "The table values are approximate and the compressible calculation is not identical to the incompressible reference.")
    println(io, "This file is a validation gate for the general-Ro convention, not a newly corrected neutral curve.")
end
println("Wrote $OUTDIR")

using Printf
using DelimitedFiles

include(joinpath(@__DIR__, "..", "src", "CompressibleBEKLinearStability.jl"))
using .CompressibleBEKLinearStability

const OUTDIR = normpath(joinpath(@__DIR__, "..", "results", "compressible_bek_reference_Mr03_20260911"))
const refs = [
    (flow="von_Karman", Ro=-1.0, alpha=0.38482, beta=0.07759, R=285.36),
    (flow="Ekman",      Ro= 0.0, alpha=0.27750, beta=0.13200, R=116.25),
    (flow="Bodewadt",   Ro= 1.0, alpha=0.52310, beta=0.12700, R=27.26),
]

rows = NamedTuple[]
for ref in refs
    p = PhysicalParameters(Mr=0.3, Tw=1.0, Ro=ref.Ro, R=ref.R,
        beta=ref.beta, omega=0.0, gamma=1.4, Pr=0.72)
    n = NumericalParameters(N=69, domain=:rational, ymax=40.0,
        target=ComplexF64(ref.alpha), neigs=4, tol=1e-11,
        maxit=800, balance=true, regularization=0.0)
    @printf("%s: Ro=%+.1f R=%.2f beta=%.5f target alpha=%.5f\n",
        ref.flow, ref.Ro, ref.R, ref.beta, ref.alpha)
    c = prepare_case(p, n; allow_unverified_ro=(ref.Ro != -1.0))
    audit = audit_expanded_operator(c)
    result = solve_mode(c)
    j = argmin(abs.(real.(result.alpha) .- ref.alpha))
    row = (flow=ref.flow, Ro=ref.Ro, Mr=p.Mr, R=ref.R, beta=ref.beta,
        alpha_ref=ref.alpha, alpha_real=real(result.alpha[j]),
        alpha_imag=imag(result.alpha[j]),
        abs_error=abs(real(result.alpha[j])-ref.alpha),
        rel_error=abs(real(result.alpha[j])-ref.alpha)/ref.alpha,
        backward_error=result.backward_errors[j],
        operator_max=maximum(audit.relative), Ma=c.Ma,
        evidence_scope=String(c.evidence_scope))
    push!(rows, row)
    @printf("  selected alpha=%.12f%+.12fi abs=%.4e rel=%.4e berr=%.3e\n",
        row.alpha_real, row.alpha_imag, row.abs_error, row.rel_error,
        row.backward_error)
end

mkpath(OUTDIR)
open(joinpath(OUTDIR, "comparison.csv"), "w") do io
    println(io, join(keys(rows[1]), ','))
    for row in rows
        println(io, join(values(row), ','))
    end
end
open(joinpath(OUTDIR, "README.txt"), "w") do io
    println(io, "Compressible BEK comparison at Mr=0.3 against supplied incompressible reference table")
    println(io, "Tw=1, omega=0, gamma=1.4, Pr=0.72, N=69, rational map ymax=40")
    println(io, "Reference values are rounded and the Ro != -1 cases use the explicitly unverified general-Ro extension.")
end
println("Wrote $OUTDIR")

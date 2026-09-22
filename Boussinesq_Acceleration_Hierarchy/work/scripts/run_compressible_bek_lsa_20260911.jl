using Printf

include(joinpath(@__DIR__, "..", "src", "CompressibleBEKLinearStability.jl"))
using .CompressibleBEKLinearStability

function options(args)
    parsed = Dict{String,String}()
    for arg in args
        startswith(arg, "--") || throw(ArgumentError("expected --name=value, got $arg"))
        pair = split(arg[3:end], '='; limit=2)
        length(pair) == 2 || throw(ArgumentError("expected --name=value, got $arg"))
        parsed[pair[1]] = pair[2]
    end
    return parsed
end

value(opts, name, default, ::Type{T}) where {T} = haskey(opts, name) ? parse(T, opts[name]) : default

opts = options(ARGS)
p = PhysicalParameters(
    Tw=value(opts, "Tw", 1.0, Float64),
    Mr=value(opts, "Mr", 0.1, Float64),
    Ro=value(opts, "Ro", -1.0, Float64),
    R=value(opts, "R", 285.36, Float64),
    beta=value(opts, "beta", 0.07759, Float64),
    omega=value(opts, "omega", 0.0, Float64),
    gamma=value(opts, "gamma", 1.4, Float64),
    Pr=value(opts, "Pr", 0.72, Float64),
)
n = NumericalParameters(
    N=value(opts, "N", 69, Int),
    domain=Symbol(get(opts, "domain", "rational")),
    ymax=value(opts, "ymax", 40.0, Float64),
    target=ComplexF64(value(opts, "target_real", 0.38482, Float64),
                      value(opts, "target_imag", 0.0, Float64)),
    neigs=value(opts, "neigs", 2, Int),
    tol=value(opts, "tol", 1.0e-11, Float64),
    maxit=value(opts, "maxit", 800, Int),
    balance=get(opts, "balance", "true") == "true",
    regularization=value(opts, "regularization", 0.0, Float64),
)
allow_unverified_ro = get(opts, "allow_unverified_ro", "false") == "true"
out = get(opts, "out", joinpath(
    @__DIR__, "..", "results", "compressible_bek_lsa_20260911",
    @sprintf("Ro_%+.6f_Mr_%.4f_Tw_%.4f_R_%.3f_N_%d", p.Ro, p.Mr, p.Tw, p.R, n.N),
))

@printf("Preparing case: Ro=%g Tw=%g Mr=%g R=%g beta=%g omega=%g N=%d\n",
        p.Ro, p.Tw, p.Mr, p.R, p.beta, p.omega, n.N)
case = prepare_case(p, n; allow_unverified_ro=allow_unverified_ro)
@printf("Internal Ma=Mr/R=%.12g; coordinate=similarity y; scope=%s\n",
        case.Ma, String(case.evidence_scope))

audit = audit_expanded_operator(case)
@printf("Operator audit relative errors: A0=%.3e A1=%.3e A2=%.3e\n", audit.relative...)

result = solve_mode(case)
checks = validate_result(result)
for j in eachindex(result.alpha)
    @printf("mode %d: alpha=%.15g%+.15gi backward_error=%.3e\n",
            j, real(result.alpha[j]), imag(result.alpha[j]), result.backward_errors[j])
end
@printf("Base BC errors: wall(F,G,H,T)=(%.3e,%.3e,%.3e,%.3e), far(F,G-1,T-1)=(%.3e,%.3e,%.3e)\n",
        checks.base_bc.wall_F, checks.base_bc.wall_G, checks.base_bc.wall_H,
        checks.base_bc.wall_T, checks.base_bc.far_F, checks.base_bc.far_G,
        checks.base_bc.far_T)
write_result(result, out)
@printf("Wrote %s\n", abspath(out))

using Printf
include(joinpath(@__DIR__, "..", "src", "CompressibleBEKLinearStability.jl"))
using .CompressibleBEKLinearStability

out = normpath(joinpath(@__DIR__, "..", "results", "compressible_bek_general_ro_operator_audit_20260911"))
mkpath(out)
ros = [-1.0, -0.5, 0.0, 0.5, 1.0]
open(joinpath(out, "directional_qep_audit.csv"), "w") do io
    println(io, "Ro,relative_directional_error,operator_max_relative_error,backward_error")
    for Ro in ros
        p = PhysicalParameters(Mr=0.01, Tw=1.0, Ro=Ro,
                               R=285.36, beta=0.07759, omega=0.0)
        n = NumericalParameters(N=49, target=0.4 + 0.0im, neigs=1,
                                tol=1e-10, regularization=0.0)
        c = prepare_case(p, n; allow_unverified_ro=true)
        d = directional_qep_test(c; alpha=0.4 + 0.02im, h=1e-6)
        a = audit_expanded_operator(c)
        r = solve_mode(c)
        @printf(io, "%.8g,%.8e,%.8e,%.8e\n", Ro, d.relative_error,
                maximum(a.relative), r.backward_errors[1])
        @printf("Ro=%+.1f directional=%.3e expanded=%.3e backward=%.3e\n",
                Ro, d.relative_error, maximum(a.relative), r.backward_errors[1])
    end
end
println("Wrote $out")

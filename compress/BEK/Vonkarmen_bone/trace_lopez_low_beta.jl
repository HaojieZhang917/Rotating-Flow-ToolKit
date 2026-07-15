include(joinpath(@__DIR__, "LopezStability.jl"))
include(joinpath(@__DIR__, "LopezBaseflow.jl"))
include(joinpath(@__DIR__, "CRD_STA.jl"))

using .LopezStability
using .LopezBaseflow
using .CRD_BF
using BSplineKit
using Printf

Tw = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.05
beta_end = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.0301
profile = solve_lopez_baseflow(Tw)
D, D2, z_cheb = CRD_BF.Cheb(99)
sample(values) = begin
    spline = BSplineKit.interpolate(profile.z, values, BSplineOrder(4))
    vec([spline(point) for point in z_cheb])
end
F, G, H, T = sample(profile.F), sample(profile.G), sample(profile.H), sample(profile.T)

output_path = joinpath(@__DIR__, "lopez_neutral_curves", "low_beta_Tw1p05.dat")
open(output_path, "w") do stream
    println(stream, "Variables=\"R\" \"beta\" \"alpha_r\" \"alpha_i\" \"residual\" \"overlap\"")
    println(stream, "Zone T=\"Tw=$Tw low beta\"")
    result = track_neutral_curve(
        F, G, H, T, 0.0, 99, D, D2;
        beta_start=0.0465,
        beta_end=beta_end,
        beta_step=-0.0004,
        R_start=352.049318,
        alpha_start=0.1461407834 + 0im,
        initial_R_delta=-1.72,
        initial_alpha_delta=-0.00272 + 0im,
        R_step=1.0,
        neutral_tol=1.0e-7,
        num_candidates=2,
        on_point=(point, index) -> begin
            @printf(
                stream, "%.16e %.16e %.16e %.16e %.16e %.16e\n",
                point.R, point.beta, real(point.alpha), imag(point.alpha),
                point.residual, point.overlap,
            )
            flush(stream)
        end,
        verbose=true,
    )
    result.failure === nothing || error("low-beta trace failed: $(result.failure)")
end

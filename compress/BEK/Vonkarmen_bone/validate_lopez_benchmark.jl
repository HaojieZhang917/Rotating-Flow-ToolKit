include(joinpath(@__DIR__, "LopezStability.jl"))
include(joinpath(@__DIR__, "LopezBaseflow.jl"))
include(joinpath(@__DIR__, "CRD_STA.jl"))
include(joinpath(@__DIR__, "Stability.jl"))

using .LopezStability
using .LopezBaseflow
using .CRD_BF
using .CRC_STA
using Printf

const PR = 0.72
const OMEGA_BENCHMARK = 0.0
const BENCHMARKS = (
    (label="Type I", R=285.36, beta=0.07759, alpha=0.38482),
    (label="Type II", R=440.88, beta=0.04672, alpha=0.13228),
)

function sample_profile(profile, points)
    continuous_solution = if profile.raw_solution isa AbstractDict
        profile.raw_solution["sol"]
    else
        getproperty(profile.raw_solution, :sol)
    end
    state = Matrix{Float64}(continuous_solution(vec(points)))
    return vec(state[3, :]), vec(state[5, :]), vec(state[1, :]), vec(state[7, :])
end

function baseflow_residuals(F, G, H, T, z, D, D2)
    chi = 2 .- T
    Fp = D * F
    Gp = D * G
    Hp = D * H
    Tp = D * T
    interior = findall(point -> 1.0e-8 < point < 20.0, z)
    return (
        F=maximum(abs, (D2 * F .- chi .* (F.^2 .+ H .* Fp .- (1 .- G).^2))[interior]),
        G=maximum(abs, (D2 * G .- chi .* (2 .* F .* G .+ H .* Gp .- 2 .* F))[interior]),
        H=maximum(abs, (Hp .+ 2 .* F)[interior]),
        T=maximum(abs, (D2 * T .- PR .* H .* Tp)[interior]),
    )
end

profile = solve_lopez_baseflow(1.0; tol=2.0e-10, nout=4001)
grid_name = isempty(ARGS) ? "malik-z20" : only(ARGS)
grid_builder = if grid_name == "current-z40"
    CRD_BF.Cheb
elseif grid_name == "malik-z20"
    CRC_STA.cheb_points
else
    error("grid must be current-z40 or malik-z20")
end
println("Lopez isothermal rotating-disk benchmarks from Malik (1986)")
println("grid=$grid_name")
for benchmark in BENCHMARKS
    println(
        "$(benchmark.label): R=$(benchmark.R) beta=$(benchmark.beta) " *
        "alpha_ref=$(benchmark.alpha) omega_D=$OMEGA_BENCHMARK",
    )
end

finest_results = Dict{String, NamedTuple}()
for N in (39, 49, 59, 69, 99)
    D, D2, z = grid_builder(N)
    F, G, H, T = sample_profile(profile, z)
    residuals = baseflow_residuals(F, G, H, T, z, D, D2)
    @printf(
        "N=%3d base_residuals_z_lt_20=(%.2e,%.2e,%.2e,%.2e)\n",
        N, residuals.F, residuals.G, residuals.H, residuals.T,
    )
    for benchmark in BENCHMARKS
        values, _, info = eigsol_lopez(
            F, G, H, T,
            benchmark.R, OMEGA_BENCHMARK, benchmark.beta,
            N, D, D2, complex(benchmark.alpha), 1;
            Pr=PR, frame=:disk, tol=1.0e-13, maxit=1200,
            return_info=true,
        )
        alpha = values[1]
        relative_real_error = abs(real(alpha) - benchmark.alpha) / benchmark.alpha
        @printf(
            "  %-7s alpha=(%.12f,%+.12e) real_error=%+.3e (%.3f%%) pencil_residual=%.3e\n",
            benchmark.label, real(alpha), imag(alpha),
            real(alpha) - benchmark.alpha, 100relative_real_error,
            info.residuals[1],
        )
        if N == 99
            finest_results[benchmark.label] = (
                alpha=alpha,
                relative_real_error=relative_real_error,
                residual=info.residuals[1],
            )
        end
    end
end

println("Acceptance checks at N=99")
for benchmark in BENCHMARKS
    result = finest_results[benchmark.label]
    passed = (
        result.relative_real_error < 5.0e-3 &&
        abs(imag(result.alpha)) < 1.0e-3 &&
        result.residual < 1.0e-10
    )
    println(
        "  $(benchmark.label): $(passed ? "PASS" : "FAIL") " *
        "(real error < 0.5%, |alpha_i| < 1e-3, residual < 1e-10)",
    )
    passed || error("$(benchmark.label) benchmark did not pass")
end

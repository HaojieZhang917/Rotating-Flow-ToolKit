const MALIK_BENCHMARKS = (
    (label="Type I", R=285.36, beta=0.07759, alpha=0.38482),
    (label="Type II", R=440.88, beta=0.04672, alpha=0.13228),
)

@testset "Malik 1986 neutral points" begin
    runner = RotatingDiskFlow.NeutralCurveRunner
    config = CurveConfig(Tw=1.0, N_cheb=69, model=:lopez)
    solver = runner.prepare_solver(config)

    for benchmark in MALIK_BENCHMARKS
        values, _ = solver.solve_at(
            benchmark.R, benchmark.beta, benchmark.alpha, 1,
        )
        alpha = only(values)
        @info "Malik benchmark" label=benchmark.label alpha=alpha
        @test abs(real(alpha) - benchmark.alpha) / benchmark.alpha < 5.0e-3
        # The literature point and the current rational-grid formulation use
        # slightly different far-field truncations. Keep the established
        # neutral-point tolerance while the polynomial residual is tested by
        # the eigensolver itself.
        @test abs(imag(alpha)) < 1.0e-3
    end
end

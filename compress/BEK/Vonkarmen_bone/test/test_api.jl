@testset "Package API and curve files" begin
    @test isdefined(RotatingDiskFlow, :NeutralCurveRunner)
    @test isdefined(RotatingDiskFlow, :LopezBaseflow)
    @test isdefined(RotatingDiskFlow, :LopezStability)
    @test isdefined(RotatingDiskFlow, :SutherlandMarching)

    runner = RotatingDiskFlow.NeutralCurveRunner
    config = CurveConfig(Tw=1.10)
    @test runner.checked_config(config) === config
    @test basename(config.output_dir) == "neutral_curve_batch"
    @test basename(dirname(config.output_dir)) == "Vonkarmen_bone"
    @test occursin("model=lopez", runner.case_tag(config))

    compressible = CurveConfig(
        Tw=1.10,
        model=:compressible,
        property_perturbations=false,
        base_property_variation=false,
    )
    tag = runner.case_tag(compressible)
    @test occursin("propPert=off", tag)
    @test occursin("baseProp=frozen", tag)

    @test_throws ArgumentError runner.checked_config(CurveConfig(Tw=1.0, N_cheb=19))
    @test_throws ArgumentError runner.checked_config(CurveConfig(Tw=1.0, model=:invalid))
    @test_throws ArgumentError runner.checked_config(CurveConfig(Tw=1.0, beta_step=0.0))

    mktempdir() do directory
        local_config = CurveConfig(
            Tw=1.10,
            output_dir=directory,
            min_valid_points=8,
        )
        beta = collect(range(0.04, step=8.0e-4, length=8))
        data = hcat(
            zeros(8),
            collect(range(300.0, 307.0, length=8)),
            beta,
            fill(0.2, 8),
            zeros(8),
            fill(0.2, 8),
            zeros(8),
        )
        path = joinpath(directory, "curve.dat")
        runner.write_curve(path, local_config, data)
        result = validate_curve_file(path, local_config)
        @test result.ok
        @test result.data == data

        data[4, 3] = data[3, 3]
        @test !runner.validate_curve_data(data, local_config).ok
    end
end

@testset "Numerical utilities" begin
    continuation = RotatingDiskFlow.NeutralContinuation
    @test continuation.refined_beta_step(8.0e-4, 5.0e-5) == 4.0e-4
    @test continuation.refined_beta_step(1.0e-4, 5.0e-5) == 5.0e-5
    @test continuation.refined_beta_step(5.0e-5, 5.0e-5) === nothing
    @test continuation.refined_beta_step(8.0e-4, 5.0e-5; enabled=false) === nothing
    @test continuation.recovered_beta_step(2.0e-4, 8.0e-4, 4, 4) == 4.0e-4
    @test continuation.recovered_beta_step(2.0e-4, 8.0e-4, 3, 4) == 2.0e-4

    D, D2, z = RotatingDiskFlow.CRD_BF.Cheb(24; domain=:finite, ymax=20.0)
    points = vec(z)
    @test size(D) == (25, 25)
    @test size(D2) == (25, 25)
    @test issorted(points)
    @test maximum(abs, D * ones(25)) < 1.0e-12
    @test maximum(abs, D * points .- 1.0) < 1.0e-11
    @test maximum(abs, D2 * (points .^ 2) .- 2.0) < 1.0e-9

    marching = RotatingDiskFlow.SutherlandMarching
    params = marching.Params(verbose=false)
    @test marching.suth_mu(1.0, params) ≈ 1.0
    @test marching.suth_rho(1.0) ≈ 1.0
    @test marching.suth_mu(1.5, params) > marching.suth_mu(1.0, params)
    @test marching.suth_rho(1.5) < marching.suth_rho(1.0)
    @test marching.radial_derivative_bdf2(9.0, 4.0, 1.0, 1.0) ≈ 6.0
end

using Test
using LinearAlgebra
using Random

include(joinpath(@__DIR__, "LopezStability.jl"))
include(joinpath(@__DIR__, "Stability.jl"))

using .LopezStability
using .CRC_STA

const C = ComplexF64

function test_derivative_matrices(n)
    rng = MersenneTwister(271828)
    D = randn(rng, n, n)
    D .-= Diagonal(vec(sum(D; dims=2)))
    return D, D * D
end

function split_fields(q, n)
    return (
        q[1:n],
        q[n+1:2n],
        q[2n+1:3n],
        q[3n+1:4n],
        q[4n+1:5n],
    )
end

# Semidiscrete nonlinear Lopez residual at one Fourier phase. The basic
# velocity retains its exact local radial dependence, while each disturbance
# derivative follows exp(i(alpha*(r-R) + beta*R*theta - omega_I*t)).
function nonlinear_fourier_residual(
    epsilon, q, F, GL, H, T, Pr, D, D2, R, alpha, beta, omega_i,
)
    n = length(F)
    u, v, w, theta, p = split_fields(q, n)
    Q = 1 .- GL

    U = F .+ epsilon .* u
    V = Q .+ epsilon .* v
    W = H ./ R .+ epsilon .* w
    temperature = T .+ epsilon .* theta
    chi = 2 .- temperature

    Ur = F ./ R .+ epsilon .* (im * alpha .* u)
    Utheta = epsilon .* (im * beta .* u)
    Uz = D * F .+ epsilon .* (D * u)

    Vr = Q ./ R .+ epsilon .* (im * alpha .* v)
    Vtheta = epsilon .* (im * beta .* v)
    Vz = D * Q .+ epsilon .* (D * v)

    Wr = epsilon .* (im * alpha .* w)
    Wtheta = epsilon .* (im * beta .* w)
    Wz = (D * H) ./ R .+ epsilon .* (D * w)

    acceleration_r = U .* Ur .+ V .* Utheta .+ W .* Uz .- V.^2 ./ R
    acceleration_theta = U .* Vr .+ V .* Vtheta .+ W .* Vz .+ U .* V ./ R
    acceleration_z = U .* Wr .+ V .* Wtheta .+ W .* Wz

    k2 = alpha^2 + beta^2
    lap_u = D2 * F .+ epsilon .* (D2 * u .- k2 .* u)
    lap_v = D2 * Q .+ epsilon .* (D2 * v .- k2 .* v)
    lap_w = (D2 * H) ./ R .+ epsilon .* (D2 * w .- k2 .* w)

    residual_u = (
        -im * omega_i * epsilon .* u .+ chi .* acceleration_r .+
        epsilon .* (im * alpha .* p) .- lap_u ./ R
    )
    residual_v = (
        -im * omega_i * epsilon .* v .+ chi .* acceleration_theta .+
        epsilon .* (im * beta .* p) .- lap_v ./ R
    )
    residual_w = (
        -im * omega_i * epsilon .* w .+ chi .* acceleration_z .+
        epsilon .* (D * p) .- lap_w ./ R
    )

    temperature_z = D * T .+ epsilon .* (D * theta)
    thermal_advection = (
        U .* (epsilon .* im * alpha .* theta) .+
        V .* (epsilon .* im * beta .* theta) .+
        W .* temperature_z
    )
    lap_temperature = D2 * T .+ epsilon .* (D2 * theta .- k2 .* theta)
    residual_temperature = (
        -im * omega_i * epsilon .* theta .+ thermal_advection .-
        lap_temperature ./ (Pr * R)
    )

    residual_continuity = (
        Ur .+ U ./ R .+ epsilon .* (im * beta .* v) .+ D * W
    )
    return vcat(
        residual_u, residual_v, residual_w,
        residual_temperature, residual_continuity,
    )
end

function finite_difference_jacobian(
    F, GL, H, T, Pr, D, D2, R, alpha, beta, omega_i; step=2.0e-6,
)
    nstate = 5length(F)
    jacobian = zeros(C, nstate, nstate)
    direction = zeros(C, nstate)
    for column in 1:nstate
        fill!(direction, 0)
        direction[column] = 1
        plus = nonlinear_fourier_residual(
            step, direction, F, GL, H, T, Pr, D, D2,
            R, alpha, beta, omega_i,
        )
        minus = nonlinear_fourier_residual(
            -step, direction, F, GL, H, T, Pr, D, D2,
            R, alpha, beta, omega_i,
        )
        jacobian[:, column] .= (plus .- minus) ./ (2step)
    end
    return jacobian
end

relative_error(left, right) = norm(left - right) / max(norm(right), eps(Float64))

@testset "Lopez generalized-Boussinesq matrices" begin
    n = 6
    N = n - 1
    D, D2 = test_derivative_matrices(n)
    z = collect(range(0.0, 3.0; length=n))
    F = 0.18 .* exp.(-0.7 .* z) .- 0.025 .* exp.(-2.0 .* z)
    GL = 1 .- exp.(-0.9 .* z)
    H = -0.72 .* (1 .- exp.(-0.6 .* z))
    T = 1 .+ 0.12 .* exp.(-0.8 .* z)
    Pr = 0.72
    R = 285.36
    beta = 0.07759
    omega_i = 0.031

    L0, L1, L2, _ = lopez_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_i; frame=:inertial,
    )
    J0 = finite_difference_jacobian(
        F, GL, H, T, Pr, D, D2, R, 0.0, beta, omega_i,
    )
    Jplus = finite_difference_jacobian(
        F, GL, H, T, Pr, D, D2, R, 1.0, beta, omega_i,
    )
    Jminus = finite_difference_jacobian(
        F, GL, H, T, Pr, D, D2, R, -1.0, beta, omega_i,
    )
    L1_fd = (Jplus - Jminus) ./ 2
    L2_fd = (Jplus + Jminus) ./ 2 - J0

    @test relative_error(L0, J0) < 2.0e-9
    @test relative_error(L1, L1_fd) < 2.0e-9
    @test relative_error(L2, L2_fd) < 2.0e-9

    alpha = 0.38402 + 0.002im
    omega_disk = -0.014 + 0.003im
    L0_disk, L1_disk, L2_disk, _ = lopez_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_disk; frame=:disk,
    )
    L0_inertial, L1_inertial, L2_inertial, _ = lopez_spatial_matrices(
        F, GL, H, T, Pr, D, D2, R, beta, omega_disk + beta;
        frame=:inertial,
    )
    @test L0_disk == L0_inertial
    @test L1_disk == L1_inertial
    @test L2_disk == L2_inertial

    A, B, _ = lopez_temporal_matrices(
        F, GL, H, T, Pr, D, D2, R, alpha, beta; frame=:disk,
    )
    temporal_operator = A - omega_disk .* B
    spatial_operator = L0_disk + alpha .* L1_disk + alpha^2 .* L2_disk
    @test relative_error(temporal_operator, spatial_operator) < 5.0e-15

    T_iso = ones(n)
    old_coefficients = CRC_STA.Spatial_mode_BEK(
        F, -GL, H, T_iso, Pr, N, D, D2, R,
    )
    old_L0, old_L1, old_L2 = CRC_STA.assemble_mat(
        old_coefficients, D, D2, beta, real(omega_disk), R,
    )
    new_L0, new_L1, new_L2, _ = lopez_spatial_matrices(
        F, GL, H, T_iso, Pr, D, D2, R, beta, real(omega_disk);
        frame=:disk,
    )
    hydro = vcat(collect(1:3n), collect(4n+1:5n))
    @test relative_error(new_L0[hydro, hydro], old_L0[hydro, hydro]) < 5.0e-15
    @test relative_error(new_L1[hydro, hydro], old_L1[hydro, hydro]) < 5.0e-15
    @test relative_error(new_L2[hydro, hydro], old_L2[hydro, hydro]) < 5.0e-15

    removed = LopezStability.homogeneous_boundary_indices(n)
    @test removed == [1, n, n + 1, 2n, 2n + 1, 3n, 3n + 1, 4n, 5n]
    @test 4n + 1 ∉ removed
    @test 5n ∈ removed
end

#!/usr/bin/env julia

"""Independent algebra/frame checks for docs/free_stream_vortical_rotating_model.md.

Run: julia --project=. scripts/verify_rotating_vortical_models.jl
This verifies analytical incident fields, not a boundary-layer or cavity solver.
Only Julia standard libraries are used.
"""

using LinearAlgebra
using Test
using Printf

const EZ = [0.0, 0.0, 1.0]
const JZ = [0.0 -1.0 0.0; 1.0 0.0 0.0; 0.0 0.0 0.0]

rotation(phi) = [cos(phi) -sin(phi) 0.0;
                 sin(phi)  cos(phi) 0.0; 0.0 0.0 1.0]

function transverse_basis(k)
    kh = hypot(k[1], k[2])
    e1 = kh > 1e-12 * norm(k) ? [-k[2], k[1], 0.0] / kh : [1.0, 0.0, 0.0]
    e2 = cross(k / norm(k), e1)
    return e1, e2
end

function core_mode(k; omega_core, viscosity, branch, axial_velocity=0.0)
    e1, e2 = transverse_basis(k)
    amplitude = (e1 - im * branch * e2) / sqrt(2)
    intrinsic = branch * 2omega_core * k[3] / norm(k) - im * viscosity * dot(k, k)
    frequency = k[3] * axial_velocity + intrinsic
    pressure = 2im * omega_core * dot(k, cross(EZ, amplitude)) / dot(k, k)
    return (; amplitude, frequency, pressure, intrinsic)
end

function disk_amplitudes(k0, mode, t, omega_core, omega_disk)
    Q = rotation((omega_core - omega_disk) * t)
    phase_time = exp(-im * mode.frequency * t)
    return Q * k0, (Q * mode.amplitude) * phase_time, mode.pressure * phase_time
end

function cylindrical_field(r, theta, z, t, k0, mode, omega_core, omega_disk)
    k, a, pressure = disk_amplitudes(k0, mode, t, omega_core, omega_disk)
    position = [r * cos(theta), r * sin(theta), z]
    phase = exp(im * dot(k, position))
    er = [cos(theta), sin(theta), 0.0]
    et = [-sin(theta), cos(theta), 0.0]
    return [dot(er, a), dot(et, a), a[3], pressure] * phase
end

function besselj_integer(n, x)
    term = (x / 2)^n / factorial(n)
    value = term
    for j in 1:160
        term *= -(x^2 / 4) / (j * (n + j))
        value += term
        abs(term) <= eps(Float64) * max(abs(value), eps(Float64)) && return value
    end
    error("Bessel series did not converge")
end

function bessel_field(r, theta, z, t; m, kh, kz, omega_core, omega_disk, W, branch)
    sigma = branch * 2omega_core * kz / hypot(kh, kz)
    frequency = m * (omega_core - omega_disk) + kz * W + sigma
    pressure = besselj_integer(abs(m), kh * r)
    pressure_r = abs(m) / r * pressure - kh * besselj_integer(abs(m) + 1, kh * r)
    denominator = 4omega_core^2 - sigma^2
    u = im * (sigma * pressure_r - 2omega_core * m / r * pressure) / denominator
    v = (2omega_core * pressure_r - sigma * m / r * pressure) / denominator
    w = kz * pressure / sigma
    return [u, v, w, pressure] * exp(im * (m * theta + kz * z - frequency * t))
end

function cylindrical_residual_fd(field, x; omega_core, omega_disk, viscosity, W)
    q = field(x)
    h = 2e-4
    derivatives = Vector{ComplexF64}[]
    second_derivatives = Vector{ComplexF64}[]
    for j in 1:4
        shift = zeros(4)
        shift[j] = h
        plus, minus = field(x + shift), field(x - shift)
        push!(derivatives, (plus - minus) / (2h))
        push!(second_derivatives, (plus - 2q + minus) / h^2)
    end
    dr, dtheta, dz, dt = derivatives
    drr, dthetatheta, dzz = second_derivatives[1:3]
    r = x[1]
    lap = drr + dr / r + dthetatheta / r^2 + dzz
    transport = dt + (omega_core - omega_disk) * dtheta + W * dz
    radial = transport[1] - 2omega_core * q[2] + dr[4] -
             viscosity * (lap[1] - q[1] / r^2 - 2dtheta[2] / r^2)
    azimuthal = transport[2] + 2omega_core * q[1] + dtheta[4] / r -
                viscosity * (lap[2] - q[2] / r^2 + 2dtheta[1] / r^2)
    axial = transport[3] + dz[4] - viscosity * lap[3]
    continuity = dr[1] + q[1] / r + dtheta[2] / r + dz[3]
    return norm([radial, azimuthal, axial, continuity])
end

function main()
    BLAS.set_num_threads(1)
    diagnostics = Pair{String, Float64}[]
    max_momentum = max_divergence = max_frame = max_k_transport = 0.0
    @testset "Rotating vortical/inertial incident-model checks" begin
        for omega_core in (0.0, 0.37, 1.0), omega_disk in (0.0, 0.8, 1.0, 1.4),
            branch in (-1, 1), k0 in ([0.8, -0.5, 1.1], [0.6, 0.9, 0.0], [0.0, 0.0, 1.2])
            viscosity, W = 0.013, -0.17
            mode = core_mode(k0; omega_core, viscosity, branch, axial_velocity=W)
            k2 = dot(k0, k0)
            a, pi = mode.amplitude, mode.pressure
            momentum = (-im * (mode.frequency - k0[3] * W) + viscosity * k2) * a +
                       2omega_core * cross(EZ, a) + im * k0 * pi
            divergence = abs(dot(k0, a))
            max_momentum = max(max_momentum, norm(momentum))
            max_divergence = max(max_divergence, divergence)
            @test norm(momentum) < 2e-12
            @test divergence < 2e-12
            @test abs(norm(a) - 1) < 2e-12
            @test abs(imag(mode.frequency) + viscosity * k2) < 2e-12

            t, dt = 0.31, 1e-5
            k, amplitude, pressure = disk_amplitudes(k0, mode, t, omega_core, omega_disk)
            kp, ap, _ = disk_amplitudes(k0, mode, t + dt, omega_core, omega_disk)
            km, am, _ = disk_amplitudes(k0, mode, t - dt, omega_core, omega_disk)
            kdot_fd, adot_fd = (kp - km) / (2dt), (ap - am) / (2dt)
            A = (omega_core - omega_disk) * JZ
            # Kelvin amplitude equations in a frame rotating at omega_disk.
            k_error = norm(kdot_fd + transpose(A) * k)
            adot_rhs = -A * amplitude - 2omega_disk * cross(EZ, amplitude) -
                       im * k * pressure - viscosity * k2 * amplitude - im * k[3] * W * amplitude
            frame_error = norm(adot_fd - adot_rhs)
            max_frame = max(max_frame, frame_error)
            max_k_transport = max(max_k_transport, k_error)
            @test k_error < 2e-8
            @test frame_error < 2e-8
            if omega_core == 0
                @test abs(real(mode.intrinsic)) < 2e-12
                @test abs(pressure) < 2e-12
            end
        end

        # Pure cylindrical harmonic: Doppler sign checked by the phase mapping.
        for m in (-4, 0, 3), omega_disk in (0.0, 0.9, 1.3)
            theta, t, omega_lab = 0.47, 0.29, 0.72
            omega_disk_frame = omega_lab - m * omega_disk
            inertial_value = exp(im * (m * (theta + omega_disk * t) - omega_lab * t))
            rotating_value = exp(im * (m * theta - omega_disk_frame * t))
            @test abs(inertial_value - rotating_value) < 2e-12
        end

        # General affine shear: differentiate k·a=0 using the derived pressure.
        A = [0.12 -0.44 0.07; 0.31 -0.08 0.19; -0.05 0.02 -0.04]
        k = [0.7, -0.2, 1.1]
        P = Matrix{Float64}(I, 3, 3) - k * transpose(k) / dot(k, k)
        a = P * ComplexF64[0.2 + 0.3im, -0.4 + 0.1im, 0.5 - 0.2im]
        omega_disk, viscosity = 0.83, 0.017
        pi = 2im * dot(k, A * a + omega_disk * cross(EZ, a)) / dot(k, k)
        kdot = -transpose(A) * k
        adot = -A * a - 2omega_disk * cross(EZ, a) - im * k * pi - viscosity * dot(k, k) * a
        constraint_rate = abs(dot(kdot, a) + dot(k, adot))
        @test constraint_rate < 2e-12
        push!(diagnostics, "affine_shear_divergence_constraint_rate" => constraint_rate)

        # Independent finite differences of the full cylindrical core equations.
        omega_core, omega_disk, viscosity, W = 0.41, 1.27, 0.021, -0.13
        k0 = [0.8, -0.4, 0.7]
        mode = core_mode(k0; omega_core, viscosity, branch=1, axial_velocity=W)
        x = [1.7, 0.63, 0.9, 0.28] # r, theta, z, t
        field(y) = cylindrical_field(y..., k0, mode, omega_core, omega_disk)
        cylindrical_residual = cylindrical_residual_fd(field, x; omega_core, omega_disk, viscosity, W)
        @test cylindrical_residual < 5e-6
        push!(diagnostics, "cylindrical_core_equations_FD_residual" => cylindrical_residual)

        # Pure integer-m Bessel waves independently check the global radial model.
        max_bessel_residual = 0.0
        for m in (-3, 0, 2), branch in (-1, 1)
            field_bessel(y) = bessel_field(y...; m, kh=0.8, kz=0.7,
                                           omega_core, omega_disk, W, branch)
            residual = cylindrical_residual_fd(field_bessel, x;
                                               omega_core, omega_disk, viscosity=0.0, W)
            @test residual < 5e-6
            max_bessel_residual = max(max_bessel_residual, residual)
        end
        push!(diagnostics, "pure_m_Bessel_core_equations_max_FD_residual" => max_bessel_residual)
    end
    append!(diagnostics, ["max_plane_wave_momentum_residual" => max_momentum,
                          "max_plane_wave_divergence" => max_divergence,
                          "max_rotating_frame_amplitude_FD_error" => max_frame,
                          "max_rotating_wavevector_FD_error" => max_k_transport])
    output = normpath(joinpath(@__DIR__, "..", "benchmark_results", "rotating_vortical_models"))
    mkpath(output)
    open(joinpath(output, "verification.txt"), "w") do io
        println(io, "Analytical incident-model verification; Julia ", VERSION)
        println(io, "These checks do not solve the wall/roughness receptivity problem.")
        for (name, value) in diagnostics
            @printf(io, "%s = %.12e\n", name, value)
        end
    end
    for (name, value) in diagnostics
        @printf("%s = %.4e\n", name, value)
    end
    println("Saved ", joinpath(output, "verification.txt"))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

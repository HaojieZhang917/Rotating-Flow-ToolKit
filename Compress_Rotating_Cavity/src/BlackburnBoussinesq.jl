module BlackburnBoussinesq

using LinearAlgebra
using NPZ
using Statistics

export BlackburnConfig, BlackburnProfile, isothermal_profile,
       fixed_pressure_profile, fixed_temperature_profile,
       pressure_tangent, trace_principal_fold

Base.@kwdef struct BlackburnConfig
    re_h::Float64 = 1000.0
    pr::Float64 = 0.72
    degree::Int = 100
    tolerance::Float64 = 2e-9
end

struct BlackburnProfile
    z::Vector{Float64}
    D::Matrix{Float64}
    fields::Matrix{Float64}  # rows H, F, G, T
    pressure::Float64
    tw::Float64
    residual::Float64
    config::BlackburnConfig
end

function chebyshev_grid(degree::Int)
    degree >= 8 || throw(ArgumentError("degree must be at least 8"))
    j = collect(0:degree)
    x = cos.(pi .* j ./ degree)
    c = vcat(2.0, ones(degree - 1), 2.0) .* ((-1.0) .^ j)
    X = repeat(x, 1, degree + 1)
    dX = X .- X'
    D = (c * (1.0 ./ c)') ./ (dX + Matrix{Float64}(I, degree + 1, degree + 1))
    D .-= Diagonal(vec(sum(D; dims=2)))
    z = 0.5 .* (1.0 .- x)
    return z, -2.0 .* D
end

pack(fields) = vcat((collect(view(fields, row, :)) for row in axes(fields, 1))...)
unpack(vector, n) = permutedims(reshape(vector, n, 4))

function set_boundary!(residual, jacobian, row, column, value)
    residual[row] = value
    jacobian[row, :] .= 0.0
    jacobian[row, column] = 1.0
end

"""
Blackburn's acceleration-consistent generalized Boussinesq closure.

The linear density factor is `chi=2-T`.  The complete local acceleration is
weighted in both steady momentum equations:

    F'' = chi * (sqrt(Re_h) H F' + Re_h (F^2-G^2+Pi))
    G'' = chi * (sqrt(Re_h) H G' + 2 Re_h F G)

This is the finite-gap steady counterpart of the Blackburn canonical closure.
"""
function system(state::AbstractVector, config::BlackburnConfig,
                fixed::Symbol, value::Float64, D::AbstractMatrix)
    n = size(D, 1)
    length(state) == 4n + 1 || throw(DimensionMismatch("thermal state length"))
    fields = unpack(view(state, 1:4n), n)
    H, F, G, T = (view(fields, row, :) for row in 1:4)
    DF, DG, DT = D * F, D * G, D * T
    D2 = D * D
    sqrt_re = sqrt(config.re_h)
    re = config.re_h
    rH, rF, rG, rT = ntuple(i -> ((i - 1)n + 1):(i*n), 4)
    pi = fixed == :temperature ? state[end] : value
    tw = fixed == :pressure ? state[end] : value
    chi = 2.0 .- T
    radial_acceleration = sqrt_re .* H .* DF .+ re .* (F.^2 .- G.^2 .+ pi)
    azimuthal_acceleration = sqrt_re .* H .* DG .+ 2re .* F .* G
    residual = vcat(
        D * H + 2sqrt_re .* F,
        D2 * F - chi .* radial_acceleration,
        D2 * G - chi .* azimuthal_acceleration,
        D2 * T - config.pr * sqrt_re .* H .* DT,
        H[end],
    )
    J = zeros(4n + 1, 4n + 1)
    eye = Matrix{Float64}(I, n, n)
    J[rH, rH] .= D
    J[rH, rF] .= 2sqrt_re .* eye
    J[rF, rH] .= -sqrt_re .* Diagonal(chi .* DF)
    J[rF, rF] .= D2 - Diagonal(chi) *
        (sqrt_re .* Diagonal(H) * D + 2re .* Diagonal(F))
    J[rF, rG] .= 2re .* Diagonal(chi .* G)
    J[rF, rT] .= Diagonal(radial_acceleration)
    J[rG, rH] .= -sqrt_re .* Diagonal(chi .* DG)
    J[rG, rF] .= -2re .* Diagonal(chi .* G)
    J[rG, rG] .= D2 - Diagonal(chi) *
        (sqrt_re .* Diagonal(H) * D + 2re .* Diagonal(F))
    J[rG, rT] .= Diagonal(azimuthal_acceleration)
    J[rT, rH] .= -config.pr * sqrt_re .* Diagonal(DT)
    J[rT, rT] .= D2 - config.pr * sqrt_re .* Diagonal(H) * D
    J[end, n] = 1.0
    set_boundary!(residual, J, first(rH), first(rH), H[1])
    set_boundary!(residual, J, first(rF), first(rF), F[1])
    set_boundary!(residual, J, last(rF), last(rF), F[end])
    set_boundary!(residual, J, first(rG), first(rG), G[1] - 1)
    set_boundary!(residual, J, last(rG), last(rG), G[end])
    set_boundary!(residual, J, first(rT), first(rT), T[1] - tw)
    set_boundary!(residual, J, last(rT), last(rT), T[end] - 1)
    if fixed == :temperature
        J[rF, end] .= -re .* chi
        J[first(rF), end] = 0.0
        J[last(rF), end] = 0.0
    elseif fixed == :pressure
        J[first(rT), end] = -1.0
    else
        throw(ArgumentError("fixed must be :temperature or :pressure"))
    end
    return residual, J
end

function newton(system_fn, initial, tolerance; max_iterations=24)
    state = collect(Float64, initial)
    for _ in 1:max_iterations
        residual, J = system_fn(state)
        err = norm(residual, Inf)
        err < tolerance && return state, err
        step = -(J \ residual)
        damping = 1.0
        accepted = false
        while damping >= 2.0^-18
            trial = state .+ damping .* step
            trial_err = norm(first(system_fn(trial)), Inf)
            if isfinite(trial_err) && trial_err < err
                state = trial
                accepted = true
                break
            end
            damping *= 0.5
        end
        if !accepted
            err < 10tolerance && return state, err
            error("Blackburn Newton line search failed; residual=$err")
        end
    end
    error("Blackburn Newton did not converge")
end

function solve(config::BlackburnConfig, fields::AbstractMatrix,
               parameter::Float64, fixed::Symbol, value::Float64)
    z, D = chebyshev_grid(config.degree)
    size(fields) == (4, length(z)) || throw(DimensionMismatch("seed fields"))
    solve_system = x -> system(x, config, fixed, value, D)
    state, residual = newton(solve_system, vcat(pack(fields), parameter), config.tolerance)
    pi = fixed == :temperature ? state[end] : value
    tw = fixed == :pressure ? state[end] : value
    return BlackburnProfile(z, D, unpack(view(state, 1:4length(z)), length(z)),
                           pi, tw, residual, config)
end

function linear_sample(x, y, targets)
    out = similar(targets, Float64)
    for (i, target) in pairs(targets)
        k = clamp(searchsortedlast(x, target), 1, length(x)-1)
        a = (target - x[k]) / (x[k+1] - x[k])
        out[i] = (1-a)*y[k] + a*y[k+1]
    end
    return out
end

"""Start from the validated isothermal rotor-stator reference profile."""
function isothermal_profile(config::BlackburnConfig, reference_npz::AbstractString)
    z, D = chebyshev_grid(config.degree)
    archive = NPZ.npzread(reference_npz)
    x = vec(Float64.(archive["x"]))
    fields = zeros(4, length(z))
    fields[1, :] .= linear_sample(x, vec(Float64.(archive["w0"])), z)
    fields[2, :] .= linear_sample(x, vec(Float64.(archive["u0"])), z)
    fields[3, :] .= linear_sample(x, vec(Float64.(archive["v0"])), z)
    fields[4, :] .= 1.0
    H, F, G = (view(fields, row, :) for row in 1:3)
    estimate = ((D * D) * F - sqrt(config.re_h) .* H .* (D * F)) ./ config.re_h .-
               F.^2 .+ G.^2
    pi_guess = median(estimate[3:end-2])
    return solve(config, fields, pi_guess, :temperature, 1.0)
end

fixed_pressure_profile(pi::Real, seed::BlackburnProfile) =
    solve(seed.config, seed.fields, seed.tw, :pressure, Float64(pi))

fixed_temperature_profile(tw::Real, seed::BlackburnProfile) =
    solve(seed.config, seed.fields, seed.pressure, :temperature, Float64(tw))

function pressure_tangent(profile::BlackburnProfile)
    state = vcat(pack(profile.fields), profile.tw)
    _, J = system(state, profile.config, :pressure, profile.pressure, profile.D)
    n = length(profile.z)
    chi = 2.0 .- profile.fields[4, :]
    rhs = zeros(4n + 1)
    rhs[n+2:2n-1] .= profile.config.re_h .* chi[2:end-1]
    return (J \ rhs)[end]
end

function trace_principal_fold(isothermal::BlackburnProfile;
                              coarse_step=0.003, fine_step=0.0005,
                              pressure_tolerance=1e-10,
                              max_pressure=0.25)
    branch = BlackburnProfile[isothermal]
    previous = isothermal
    # Heating proceeds toward larger Pi for this closure.  A fold is optional:
    # the routine returns `(nothing, branch)` when the scan reaches
    # max_pressure without a sign change of dTw/dPi.
    for pi in range(isothermal.pressure, max_pressure; step=coarse_step)[2:end]
        previous = fixed_pressure_profile(pi, previous)
        push!(branch, previous)
    end
    while previous.pressure < max_pressure
        pi = min(max_pressure, previous.pressure + fine_step)
        next = fixed_pressure_profile(pi, previous)
        push!(branch, next)
        g0, g1 = pressure_tangent(previous), pressure_tangent(next)
        if g0 * g1 <= 0
            lo, hi = previous, next
            for _ in 1:40
                hi.pressure - lo.pressure <= pressure_tolerance && break
                mid_pi = 0.5 * (lo.pressure + hi.pressure)
                seed = abs(mid_pi - lo.pressure) < abs(mid_pi - hi.pressure) ? lo : hi
                mid = fixed_pressure_profile(mid_pi, seed)
                if pressure_tangent(mid) * pressure_tangent(lo) <= 0
                    hi = mid
                else
                    lo = mid
                end
            end
            fold = fixed_pressure_profile(0.5 * (lo.pressure + hi.pressure), lo)
            return fold, branch
        end
        previous = next
    end
    return nothing, branch
end

end

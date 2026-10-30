module SingleDiskRadial

using LinearAlgebra

export RadialConfig, reference_profile, inlet_profile, march, interpolate_profile,
       radial_system, profile_error

"""Steady, axisymmetric, zero-pressure-gradient single-disk boundary layer.

The marching coordinate is s=log(r/r0), eta=z*sqrt(Omega/nu), and the
velocities are (Omega*r*F, Omega*r*G, sqrt(nu*Omega)*H).  This is a
non-similar boundary-layer calculation, not a time integrator or full NS.
"""
Base.@kwdef struct RadialConfig
    degree::Int = 80
    eta_max::Float64 = 20.0
    ds::Float64 = 0.02
    s_end::Float64 = log(100.0)
    method::Symbol = :bdf2
    tolerance::Float64 = 2e-10
    max_iterations::Int = 24
    reverse_tolerance::Float64 = 1e-9
end

function validate_config(config::RadialConfig)
    config.degree >= 16 || throw(ArgumentError("degree must be at least 16"))
    isfinite(config.eta_max) && config.eta_max > 0 ||
        throw(ArgumentError("eta_max must be positive and finite"))
    isfinite(config.ds) && config.ds > 0 || throw(ArgumentError("ds must be positive and finite"))
    isfinite(config.s_end) && 0 < config.s_end < 700 ||
        throw(ArgumentError("s_end must be positive, finite, and below 700"))
    config.method in (:bdf2, :be) || throw(ArgumentError("method must be :bdf2 or :be"))
    isfinite(config.tolerance) && config.tolerance > 0 ||
        throw(ArgumentError("tolerance must be positive and finite"))
    config.max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    isfinite(config.reverse_tolerance) && config.reverse_tolerance >= 0 ||
        throw(ArgumentError("reverse_tolerance must be finite and nonnegative"))
    return nothing
end

function grid(config::RadialConfig)
    validate_config(config)
    N = config.degree
    angles = pi .* (0:N) ./ N
    x = cos.(angles)
    c = vcat(2.0, ones(N - 1), 2.0) .* (-1.0) .^ (0:N)
    dx = x .- transpose(x)
    D = (c * transpose(1.0 ./ c)) ./ (dx + I)
    D .-= Diagonal(vec(sum(D; dims=2)))
    D .*= -2 / config.eta_max
    eta = config.eta_max / 2 .* (1 .- x)

    # Clenshaw-Curtis weights, including the physical eta interval scale.
    weights = zeros(N + 1)
    v = ones(N - 1)
    theta = angles[2:N]
    if iseven(N)
        weights[1] = weights[end] = 1 / (N^2 - 1)
        for k in 1:(N ÷ 2 - 1)
            v .-= 2 .* cos.(2k .* theta) ./ (4k^2 - 1)
        end
        v .-= cos.(N .* theta) ./ (N^2 - 1)
    else
        weights[1] = weights[end] = 1 / N^2
        for k in 1:((N - 1) ÷ 2)
            v .-= 2 .* cos.(2k .* theta) ./ (4k^2 - 1)
        end
    end
    weights[2:N] .= 2 .* v ./ N
    weights .*= config.eta_max / 2
    return (; config, eta, D, D2=D * D, weights)
end

function boundary_row!(residual, J, row, value)
    residual[row] = value
    if J !== nothing
        J[row, :] .= 0
        J[row, row] = 1
    end
end

"""Residual and analytic Jacobian in state order [H; F; G].

Fs=a*F+history_f and Gs=a*G+history_g.  Set a=0 and zero histories for
the von Karman reference BVP.  H(eta_max) is deliberately NOT prescribed.
"""
function radial_system(state::AbstractVector, mesh;
                       a::Real=0.0, history_f=zeros(length(mesh.eta)),
                       history_g=zeros(length(mesh.eta)), jacobian::Bool=true)
    n = length(mesh.eta)
    length(state) == 3n || throw(DimensionMismatch("state must contain H,F,G"))
    H, F, G = (view(state, (j*n+1):((j+1)*n)) for j in 0:2)
    D, D2 = mesh.D, mesh.D2
    DF, DG = D * F, D * G
    Fs, Gs = a .* F .+ history_f, a .* G .+ history_g
    residual = vcat(D * H .+ 2 .* F .+ Fs,
                    D2 * F .- H .* DF .- F.^2 .+ G.^2 .- F .* Fs,
                    D2 * G .- H .* DG .- 2 .* F .* G .- F .* Gs)
    J = jacobian ? zeros(3n, 3n) : nothing
    if jacobian
        ih, jf, kg = 1:n, (n+1):(2n), (2n+1):(3n)
        J[ih, ih] .= D
        J[ih, jf] .= (2 + a) .* Matrix{Float64}(I, n, n)
        J[jf, ih] .= -Diagonal(DF)
        J[jf, jf] .= D2 - Diagonal(H) * D - Diagonal(2 .* F .+ Fs .+ a .* F)
        J[jf, kg] .= 2 .* Diagonal(G)
        J[kg, ih] .= -Diagonal(DG)
        J[kg, jf] .= -Diagonal(2 .* G .+ Gs)
        J[kg, kg] .= D2 - Diagonal(H) * D - Diagonal((2 + a) .* F)
    end
    boundary_row!(residual, J, 1, H[1])
    boundary_row!(residual, J, n + 1, F[1])
    boundary_row!(residual, J, 2n, F[end])
    boundary_row!(residual, J, 2n + 1, G[1] - 1)
    boundary_row!(residual, J, 3n, G[end])
    return residual, J
end

function newton(system, seed, config)
    state = copy(seed)
    for iteration in 0:config.max_iterations
        residual, J = system(state, true)
        err = norm(residual, Inf)
        isfinite(err) || error("non-finite nonlinear residual")
        err <= config.tolerance && return (; state, iterations=iteration, residual=err)
        iteration == config.max_iterations && break
        direction = -(J \ residual)
        alpha = 1.0
        accepted = false
        while alpha >= 2.0^-18
            candidate = state .+ alpha .* direction
            next_residual, _ = system(candidate, false)
            next_err = norm(next_residual, Inf)
            if isfinite(next_err) && (next_err < (1 - 1e-4*alpha)*err ||
                                      next_err <= config.tolerance)
                state = candidate
                accepted = true
                break
            end
            alpha /= 2
        end
        accepted || error("Newton line search failed (residual=$err)")
    end
    error("Newton exceeded $(config.max_iterations) iterations")
end

"""Compute the main von Karman branch on the same eta grid used for marching."""
function reference_profile(config::RadialConfig=RadialConfig())
    mesh = grid(config)
    z = mesh.eta
    F = 0.442 .* z .* exp.(-z)
    G = exp.(-z)
    H = -0.884 .* (1 .- (1 .+ z) .* exp.(-z))
    F[1] = F[end] = 0
    G[1], G[end] = 1, 0
    system = (q, j) -> radial_system(q, mesh; jacobian=j)
    result = newton(system, vcat(H, F, G), config)
    n = length(z)
    H, F, G = (copy(result.state[(j*n+1):((j+1)*n)]) for j in 0:2)
    minimum(F[2:end-1]) >= -config.reverse_tolerance ||
        error("reference Newton solve did not reach the outward-flow main branch")
    return (; mesh..., H, F, G, state=result.state,
            residual=result.residual, iterations=result.iterations)
end

"""Five reproducible inlet profiles; H is not independent inlet data.

Multiplicative bumps vanish to fourth/sixth order at the wall, preserving
the low-order wall compatibility conditions.  Even the large cases retain
outward radial velocity.  Different centers/powers avoid a single-shape test.
"""
function inlet_profile(reference; case::Symbol=:similarity)
    parameters = if case == :similarity
        (af=0.0, ag=0.0, center=2.0, power=4,
         description="discrete von Karman reference")
    elseif case == :small_shape
        (af=0.02, ag=-0.02, center=1.5, power=4,
         description="small radial increase and swirl decrease")
    elseif case == :small_opposite
        (af=-0.02, ag=0.02, center=3.0, power=6,
         description="small opposite-sign, broader inlet perturbation")
    elseif case == :large_radial
        (af=1.5, ag=-0.7, center=1.8, power=4,
         description="large radial increase and swirl reduction")
    elseif case == :large_swirl
        (af=-0.65, ag=0.8, center=2.5, power=6,
         description="large radial reduction and broader swirl increase")
    else
        throw(ArgumentError("unknown inlet case: $case"))
    end
    z = reference.eta
    bump = (z ./ parameters.center).^parameters.power .*
           exp.(parameters.power .* (1 .- z ./ parameters.center))
    F = reference.F .* (1 .+ parameters.af .* bump)
    G = reference.G .* (1 .+ parameters.ag .* bump)
    F[1] = F[end] = 0
    G[1], G[end] = 1, 0
    return (; case, F, G, parameters...)
end

function profile_error(reference, F, G)
    w = reference.weights
    ef = sqrt(sum(w .* (F .- reference.F).^2) / sum(w .* reference.F.^2))
    eg = sqrt(sum(w .* (G .- reference.G).^2) / sum(w .* reference.G.^2))
    combined = sqrt(sum(w .* ((F .- reference.F).^2 .+ (G .- reference.G).^2)) /
                    sum(w .* (reference.F.^2 .+ reference.G.^2)))
    return (; error_fg=combined, error_f=ef, error_g=eg)
end

function metric(reference, F, G, H, s, station, residual, iterations, Fs)
    errors = profile_error(reference, F, G)
    computed = all(isfinite, H)
    balance = computed ? H[end] - H[1] +
                         sum(reference.weights .* (2 .* F .+ Fs)) : NaN
    return (; station, s, r_ratio=exp(s), errors...,
            H_top=computed ? H[end] : NaN,
            wall_Fp=(reference.D * F)[1], wall_Gp=(reference.D * G)[1],
            residual, iterations, min_F=minimum(F[2:end-1]),
            continuity_balance=balance)
end

"""Implicit radial marching without clipping, FLARE, or similarity projection.

BE initializes the BDF2 history.  Only F,G are prescribed at s=0; the
corresponding inlet H is not computed and is stored as NaN.  An integral
estimate of H is used solely as the first Newton guess.  A converged reverse
flow stops the calculation with status=:reverse_flow rather than changing PDEs.
"""
function march(reference, inlet; config::RadialConfig=reference.config)
    validate_config(config)
    config.degree == reference.config.degree && config.eta_max == reference.config.eta_max ||
        throw(ArgumentError("marching and reference eta grids must match"))
    n = length(reference.eta)
    length(inlet.F) == length(inlet.G) == n || throw(DimensionMismatch("inlet grid"))
    all(isfinite, inlet.F) && all(isfinite, inlet.G) || throw(ArgumentError("non-finite inlet"))
    maximum(abs.([inlet.F[1], inlet.F[end], inlet.G[1]-1, inlet.G[end]])) < 1e-10 ||
        throw(ArgumentError("inlet violates disk/far-field velocity boundary conditions"))
    steps = ceil(Int, config.s_end / config.ds)
    ds = config.s_end / steps
    F = zeros(n, steps + 1)
    G = similar(F)
    H = fill(NaN, n, steps + 1)
    F[:, 1] .= inlet.F
    G[:, 1] .= inlet.G
    Hguess = zeros(n)
    for j in 2:n
        Hguess[j] = Hguess[j-1] -
            (reference.eta[j]-reference.eta[j-1]) * (inlet.F[j]+inlet.F[j-1])
    end
    seed = vcat(Hguess, inlet.F, inlet.G)
    metrics = [metric(reference, inlet.F, inlet.G, H[:, 1], 0.0, 0, NaN, 0, zeros(n))]
    status, message = :completed, "reached requested radius ratio"
    completed = 0
    if minimum(inlet.F[2:end-1]) < -config.reverse_tolerance
        status, message = :reverse_flow, "inlet has inward radial velocity"
    else
        for k in 1:steps
            if config.method == :be || k == 1
                a = 1 / ds
                history_f = -F[:, k] ./ ds
                history_g = -G[:, k] ./ ds
            else
                a = 1.5 / ds
                history_f = (-2 .* F[:, k] .+ 0.5 .* F[:, k-1]) ./ ds
                history_g = (-2 .* G[:, k] .+ 0.5 .* G[:, k-1]) ./ ds
            end
            system = (q, j) -> radial_system(q, reference; a, history_f, history_g, jacobian=j)
            result = try
                newton(system, seed, config)
            catch err
                err isa InterruptException && rethrow()
                status = :nonlinear_failure
                message = "station $k: " * sprint(showerror, err)
                break
            end
            H[:, k+1] .= result.state[1:n]
            F[:, k+1] .= result.state[(n+1):(2n)]
            G[:, k+1] .= result.state[(2n+1):(3n)]
            Fs = a .* F[:, k+1] .+ history_f
            push!(metrics, metric(reference, F[:, k+1], G[:, k+1], H[:, k+1],
                                  k * ds, k, result.residual, result.iterations, Fs))
            completed = k
            seed = result.state
            if metrics[end].min_F < -config.reverse_tolerance
                status = :reverse_flow
                message = "converged inward radial velocity at radius ratio $(exp(k*ds)); no FLARE used"
                break
            end
        end
    end
    return (; config, reference, inlet, ds, status, message, metrics,
            s=collect(0:completed) .* ds,
            F=F[:, 1:(completed+1)], G=G[:, 1:(completed+1)], H=H[:, 1:(completed+1)])
end

"""Barycentric interpolation from a Chebyshev-Lobatto grid to physical eta points."""
function interpolate_profile(eta::AbstractVector, values::AbstractVector, points)
    n = length(eta)
    length(values) == n || throw(DimensionMismatch("profile grid"))
    bary = (-1.0) .^ (0:(n-1))
    bary[1] *= 0.5
    bary[end] *= 0.5
    output = zeros(length(points))
    for (k, point) in enumerate(points)
        first(eta) <= point <= last(eta) || throw(ArgumentError("interpolation outside eta domain"))
        nearest = argmin(abs.(eta .- point))
        if abs(eta[nearest] - point) <= 16eps(max(1.0, abs(point)))
            output[k] = values[nearest]
        else
            terms = bary ./ (point .- eta)
            output[k] = dot(terms, values) / sum(terms)
        end
    end
    return output
end

end

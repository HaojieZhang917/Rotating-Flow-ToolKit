using LinearAlgebra
using Printf
using DelimitedFiles

Base.@kwdef struct MarchParams
    Neta::Int = 64
    Re_s::Float64 = 2000.0
    R_ini::Float64 = 0.05
    R_end::Float64 = 0.25
    Nx::Int = 40
    use_bdf2::Bool = false

    # Cw is an imposed increment of the gap-averaged axial velocity H
    # relative to the no-throughflow Batchelor solution. H is scaled with
    # sqrt(Omega*nu), following Corral & Romera's convention.
    Cw::Float64 = 1e-4

    output_file::String = joinpath(@__DIR__, "throughflow_marching_profile.csv")

    newton_tol::Float64 = 1e-10
    newton_maxit::Int = 40
    fd_eps::Float64 = sqrt(eps(Float64))
    min_linesearch_alpha::Float64 = 1e-6

    # Strong negative radial velocity makes outward marching questionable.
    # Natural Batchelor flow already has a mild inward stator-side return flow,
    # so the default abort threshold is deliberately conservative.
    backflow_abort::Float64 = -0.30
    verbose::Bool = true
end

struct MarchResult
    eta::Vector{Float64}
    R::Vector{Float64}
    x::Vector{Float64}
    F::Matrix{Float64}
    G::Matrix{Float64}
    H::Matrix{Float64}
    B::Vector{Float64}
    residual::Vector{Float64}
    target_mean_H::Float64
    stopped::Bool
    stop_reason::String
end

function clenshaw_curtis_weights(N::Int)
    w = zeros(N + 1)
    if N == 0
        w[1] = 2.0
        return w
    end

    for j in 0:N
        s = 0.0
        for k in 1:floor(Int, N / 2)
            term = 2.0 / (1.0 - (2k)^2) * cos(2k * j * pi / N)
            s += (2k == N) ? 0.5 * term : term
        end
        cj = (j == 0 || j == N) ? 1.0 : 2.0
        w[j + 1] = cj / N * (1.0 + s)
    end
    return w
end

function cheb_lobatto(N::Int)
    N < 2 && error("Neta must be at least 2")

    x = [cos(pi * j / N) for j in 0:N]
    c = [j == 0 || j == N ? 2.0 : 1.0 for j in 0:N]
    c .*= [(-1.0)^j for j in 0:N]

    D = zeros(N + 1, N + 1)
    for i in 1:N + 1
        for j in 1:N + 1
            if i != j
                D[i, j] = (c[i] / c[j]) / (x[i] - x[j])
            end
        end
    end
    for i in 1:N + 1
        D[i, i] = -sum(D[i, :])
    end

    # Reverse from x in [1,-1] to eta in [0,1].
    perm = collect(N + 1:-1:1)
    x = x[perm]
    D = D[perm, perm]
    eta = (x .+ 1.0) ./ 2.0

    Deta = 2.0 .* D
    D2eta = Deta * Deta
    weights = clenshaw_curtis_weights(N)[perm] ./ 2.0

    err = norm(Deta * eta .- 1.0, Inf)
    err > 1e-8 && error("Chebyshev derivative check failed: ||D*eta-1|| = $err")
    abs(sum(weights) - 1.0) > 1e-12 && error("Quadrature weights do not sum to one")
    return eta, Deta, D2eta, weights
end

pack(F, G, H, B) = vcat(F, G, H, [B])

function unpack(q::Vector{Float64}, n::Int)
    F = q[1:n]
    G = q[n + 1:2n]
    H = q[2n + 1:3n]
    B = q[3n + 1]
    return F, G, H, B
end

function finite_difference_jacobian(f, x::Vector{Float64}, eps_scale::Float64)
    fx = f(x)
    J = zeros(length(fx), length(x))
    for j in eachindex(x)
        h = eps_scale * max(1.0, abs(x[j]))
        xp = copy(x)
        xp[j] += h
        J[:, j] = (f(xp) - fx) ./ h
    end
    return J, fx
end

function newton_solve(f, x0::Vector{Float64}, p::MarchParams; label::String = "")
    x = copy(x0)
    last_res = Inf

    for it in 1:p.newton_maxit
        J, R = finite_difference_jacobian(f, x, p.fd_eps)
        res = norm(R, Inf)
        last_res = res
        p.verbose && @printf("%s Newton %2d: ||R||_inf = %.4e\n", label, it, res)
        res < p.newton_tol && return x, res, it, true, "converged"

        dx = try
            -(J \ R)
        catch err
            return x, res, it, false, "linear solve failed: $(err)"
        end

        alpha = 1.0
        accepted = false
        while alpha >= p.min_linesearch_alpha
            trial = x .+ alpha .* dx
            trial_res = norm(f(trial), Inf)
            if isfinite(trial_res) && trial_res < res
                x = trial
                accepted = true
                break
            end
            alpha *= 0.5
        end
        accepted || return x, res, it, false, "line search failed"
    end

    return x, last_res, p.newton_maxit, false, "maximum Newton iterations reached"
end

function similarity_residual(y::Vector{Float64}, D, Re_s::Float64)
    n = size(D, 1)
    H = y[1:n]
    F = y[n + 1:2n]
    Feta = y[2n + 1:3n]
    Fetaeta = y[3n + 1:4n]
    G = y[4n + 1:5n]
    Geta = y[5n + 1:6n]

    sqrt_Re = sqrt(Re_s)
    Heta = -2.0 .* sqrt_Re .* F

    R = zeros(6n)
    R[1:n] = D * H .- Heta
    R[n + 1:2n] = D * F .- Feta
    R[2n + 1:3n] = D * Feta .- Fetaeta
    R[3n + 1:4n] = D * Fetaeta .-
        Re_s .* ((H .* Fetaeta .+ Feta .* Heta) ./ sqrt_Re .- 2.0 .* G .* Geta .+ 2.0 .* F .* Feta)
    R[4n + 1:5n] = D * G .- Geta
    R[5n + 1:6n] = D * Geta .-
        Re_s .* ((H .* Geta) ./ sqrt_Re .+ 2.0 .* F .* G)

    R[1] = H[1]
    R[n] = H[end]
    R[n + 1] = F[1]
    R[2n] = F[end]
    R[4n + 1] = G[1] - 1.0
    R[5n] = G[end]
    return R
end

function solve_similarity(eta, D, D2, weights, p::MarchParams)
    n = length(eta)
    H = zeros(n)
    F = zeros(n)
    Feta = zeros(n)
    Fetaeta = zeros(n)
    G = 1.0 .- eta
    Geta = -ones(n)
    y0 = vcat(H, F, Feta, Fetaeta, G, Geta)

    f = y -> similarity_residual(y, D, p.Re_s)
    y, res, _, ok, reason = newton_solve(f, y0, p; label = "similarity")
    ok || error("Similarity solve failed: $reason; residual = $res")

    H = y[1:n]
    F = y[n + 1:2n]
    Feta = y[2n + 1:3n]
    Fetaeta = y[3n + 1:4n]
    G = y[4n + 1:5n]

    eps_s = inv(sqrt(p.Re_s))
    B_values = eps_s^2 .* Fetaeta .- eps_s .* H .* Feta .+ G .^ 2 .- F .^ 2
    B = sum(weights .* B_values)
    return pack(F, G, H, B), res
end

function solve_discrete_stationary(q_guess, D, D2, weights, p::MarchParams)
    f = q -> marching_residual(q, q, nothing, 1, 1.0, D, D2, weights, 0.0, p)
    q, res, _, ok, reason = newton_solve(f, q_guess, p; label = "stationary")
    ok || error("Discrete stationary projection failed: $reason; residual = $res")
    return q, res
end

function inlet_perturbation(eta, D, q_base, target_mean_H, weights, p::MarchParams)
    n = length(eta)
    F0, G0, H0, B0 = unpack(q_base, n)
    current_mean = dot(weights, H0)
    dmean = target_mean_H - current_mean

    # Smooth wall-compatible shape with zero value and zero slope at both walls.
    shape = 30.0 .* eta .^ 2 .* (1.0 .- eta) .^ 2
    shape ./= dot(weights, shape)

    H = H0 .+ dmean .* shape
    eps_s = inv(sqrt(p.Re_s))
    F = F0 .- 0.5 .* eps_s .* dmean .* (D * shape)
    F[1] = 0.0
    F[end] = 0.0
    G = copy(G0)
    return pack(F, G, H, B0)
end

function marching_residual(q, q_prev, q_prevprev, step::Int, dx::Float64,
                           D, D2, weights, target_mean_H::Float64, p::MarchParams)
    n = size(D, 1)
    F, G, H, B = unpack(q, n)
    Fp, Gp, _, _ = unpack(q_prev, n)

    if !p.use_bdf2 || step == 1 || q_prevprev === nothing
        dFdx = (F .- Fp) ./ dx
        dGdx = (G .- Gp) ./ dx
    else
        Fpp, Gpp, _, _ = unpack(q_prevprev, n)
        dFdx = (3.0 .* F .- 4.0 .* Fp .+ Fpp) ./ (2.0 * dx)
        dGdx = (3.0 .* G .- 4.0 .* Gp .+ Gpp) ./ (2.0 * dx)
    end

    eps_s = inv(sqrt(p.Re_s))
    Feta = D * F
    Geta = D * G
    Heta = D * H
    Fetaeta = D2 * F
    Getaeta = D2 * G

    Rc = dFdx .+ 2.0 .* F .+ eps_s .* Heta
    Rf = F .* (F .+ dFdx) .+ eps_s .* H .* Feta .- G .^ 2 .+ B .- eps_s^2 .* Fetaeta
    Rg = F .* (G .+ dGdx) .+ eps_s .* H .* Geta .+ F .* G .- eps_s^2 .* Getaeta

    # Tau rows: solid rotor/stator walls. The inlet flux is used only to
    # construct q at R_ini. Downstream, B is determined by the second
    # impermeability condition H(eta=1)=0.
    Rc[1] = H[1]
    Rf[1] = F[1]
    Rf[end] = F[end]
    Rg[1] = G[1] - 1.0
    Rg[end] = G[end]

    Rwall = H[end]
    return vcat(Rc, Rf, Rg, [Rwall])
end

function write_output(result::MarchResult, p::MarchParams)
    rows = Float64[]
    nR = length(result.R)
    nEta = length(result.eta)

    open(p.output_file, "w") do io
        println(io, "iR,R,x,eta,F,G,H,B,meanH,minF,maxF,residual")
        for i in 1:nR
            meanH = sum(result.H[:, i]) / nEta
            minF = minimum(result.F[:, i])
            maxF = maximum(result.F[:, i])
            for j in 1:nEta
                @printf(io, "%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                        i, result.R[i], result.x[i], result.eta[j],
                        result.F[j, i], result.G[j, i], result.H[j, i],
                        result.B[i], meanH, minF, maxF, result.residual[i])
            end
        end
    end
end

function solve_marching(p::MarchParams = MarchParams())
    p.R_ini <= 0.0 && error("R_ini must be positive")
    p.R_end <= p.R_ini && error("R_end must be larger than R_ini")
    p.Nx < 1 && error("Nx must be positive")

    eta, D, D2, weights = cheb_lobatto(p.Neta)
    x_end = log(p.R_end / p.R_ini)
    dx = x_end / p.Nx
    xgrid = collect(0:p.Nx) .* dx
    Rgrid = p.R_ini .* exp.(xgrid)

    q_base, similarity_res = solve_similarity(eta, D, D2, weights, p)
    q_base, base_res = solve_discrete_stationary(q_base, D, D2, weights, p)
    F_base, G_base, H_base, B_base = unpack(q_base, length(eta))
    base_mean_H = dot(weights, H_base)
    target_mean_H = base_mean_H + p.Cw

    p.verbose && begin
        @printf("\nBase Batchelor solution:\n")
        @printf("  Re_s = %.6g, eps = %.6e\n", p.Re_s, inv(sqrt(p.Re_s)))
        @printf("  B = %.8e\n", B_base)
        @printf("  similarity residual = %.4e, stationary residual = %.4e\n", similarity_res, base_res)
        @printf("  mean(H)_base = %.8e\n", base_mean_H)
        @printf("  target mean(H) = %.8e  (Cw increment = %.8e)\n\n", target_mean_H, p.Cw)
    end

    n = length(eta)
    Fhist = zeros(n, p.Nx + 1)
    Ghist = zeros(n, p.Nx + 1)
    Hhist = zeros(n, p.Nx + 1)
    Bhist = zeros(p.Nx + 1)
    residuals = fill(NaN, p.Nx + 1)

    q0 = inlet_perturbation(eta, D, q_base, target_mean_H, weights, p)
    F0, G0, H0, B0 = unpack(q0, n)
    Fhist[:, 1] = F0
    Ghist[:, 1] = G0
    Hhist[:, 1] = H0
    Bhist[1] = B0
    residuals[1] = base_res

    q_prevprev = nothing
    q_prev = q0
    stopped = false
    stop_reason = ""
    last_index = p.Nx + 1

    for step in 1:p.Nx
        label = @sprintf("x-step %03d/%03d", step, p.Nx)
        f = q -> marching_residual(q, q_prev, q_prevprev, step, dx, D, D2, weights,
                                   target_mean_H, p)

        q_guess = copy(q_prev)
        q_new, res, _, ok, reason = newton_solve(f, q_guess, p; label)
        F, G, H, B = unpack(q_new, n)

        if !ok
            stopped = true
            stop_reason = "Newton failed at step $step, R=$(Rgrid[step + 1]): $reason; residual=$res. The last converged station was saved; the parabolized outward march is likely incompatible with the return-flow/backflow region at this radius."
            last_index = step
            break
        end

        if minimum(F) < p.backflow_abort
            stopped = true
            stop_reason = "outward marching stopped at step $step, R=$(Rgrid[step + 1]): min(F)=$(minimum(F)) < backflow_abort=$(p.backflow_abort). Strong inward radial flow makes the parabolized outward march unreliable."
            last_index = step
            break
        end

        Fhist[:, step + 1] = F
        Ghist[:, step + 1] = G
        Hhist[:, step + 1] = H
        Bhist[step + 1] = B
        residuals[step + 1] = res

        q_prevprev = q_prev
        q_prev = q_new
    end

    if stopped
        Rgrid = Rgrid[1:last_index]
        xgrid = xgrid[1:last_index]
        Fhist = Fhist[:, 1:last_index]
        Ghist = Ghist[:, 1:last_index]
        Hhist = Hhist[:, 1:last_index]
        Bhist = Bhist[1:last_index]
        residuals = residuals[1:last_index]
    end

    result = MarchResult(eta, Rgrid, xgrid, Fhist, Ghist, Hhist, Bhist, residuals,
                         target_mean_H, stopped, stop_reason)
    write_output(result, p)

    @printf("\nMarching finished.\n")
    @printf("  Saved: %s\n", p.output_file)
    @printf("  Stations saved: %d\n", length(result.R))
    @printf("  Final R = %.8e\n", result.R[end])
    @printf("  Final residual = %.4e\n", result.residual[end])
    @printf("  Final B = %.8e\n", result.B[end])
    @printf("  Final min(F) = %.8e, max(F) = %.8e\n", minimum(result.F[:, end]), maximum(result.F[:, end]))
    @printf("  Final weighted mean(H) = %.8e, target = %.8e\n",
            dot(weights, result.H[:, end]), result.target_mean_H)
    result.stopped && @printf("  Stop reason: %s\n", result.stop_reason)

    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    solve_marching()
end

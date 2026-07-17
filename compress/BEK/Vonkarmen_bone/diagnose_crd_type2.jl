using BSplineKit
using LinearAlgebra
using NonlinearEigenproblems

include("CRD_STA.jl")

const SIGMA = 0.72
const GAMMA = 1.4
const RO = -1.0
const CO = 2.0

function prepare_baseflow(; Tw=1.0, Mr=0.3)
    u0, v0, w0, f, q, _, _, _ = baseflow_var(10, RO, CO)
    H, T = T_ca(Mr, f, q, w0, GAMMA, Tw)
    return u0, v0, H, T
end

function finite_chebyshev(N, ymax)
    theta = range(0.0, pi; length=N + 1)
    eta = reshape(-cos.(theta), N + 1, 1)
    weights = [2; ones(N - 1); 2] .* (-1.0) .^ (0:N)
    eta_matrix = repeat(eta, 1, N + 1)
    delta = eta_matrix - transpose(eta_matrix)
    Deta = (weights * transpose(1 ./ weights)) ./ (delta + I)
    Deta -= diagm(vec(sum(Deta; dims=2)))
    D = (2 / ymax) .* Deta
    D2 = D^2
    y = (ymax / 2) .* (eta .+ 1)
    return D, D2, y
end

function matrices_at_resolution(base, N; coordinate="sim", R=440.88,
                                Mr=0.3, beta=0.04672, omega=0.0,
                                regularization=1e-8, grid="rational",
                                ymax=20.0)
    u0, v0, H0, T0 = base
    D, D2, y = grid == "finite" ? finite_chebyshev(N, ymax) : CRD_BF.Cheb(N)
    F, G, H, T, rho, _ = interp(u0, v0, H0, T0, y, N, coordinate)
    lambda = -(2 / 3) .* T
    kappa = T ./ SIGMA
    Ma = Mr / R

    coefficients = Spatial_mode_BEK(
        F, G, H, rho, lambda, kappa, T, SIGMA, GAMMA,
        R, Ma, N, RO, CO, D, D2,
    )
    A0, A1, A2 = assemble_mat(coefficients, D, D2, beta, omega)

    n = N + 1
    for i in 1:n
        A2[i, i] += regularization - 1e-8
    end
    A0, A1, A2 = boudary_condition(A0, A1, A2, N)

    clipped = grid == "finite" ? 0 : count(==(40.0), vec(y))
    return A0, A1, A2, clipped
end

function polynomial_residual(A0, A1, A2, alpha, vector)
    numerator = norm((A0 + alpha .* A1 + alpha^2 .* A2) * vector)
    denominator = (
        norm(A0) + abs(alpha) * norm(A1) + abs2(alpha) * norm(A2)
    ) * norm(vector)
    return numerator / max(denominator, eps(Float64))
end

function closest_finite(values, target; count=4)
    valid = findall(z -> isfinite(real(z)) && isfinite(imag(z)) && abs(z) < 5, values)
    order = sort(valid; by=i -> abs(values[i] - target))
    return first(order, min(count, length(order)))
end

function dense_companion(A0, A1, A2, target; equilibrate=false)
    m = size(A0, 1)
    column_scale = ones(Float64, m)
    if equilibrate
        row_measure = vec(maximum(abs.(hcat(A0, A1, A2)); dims=2))
        row_scale = 1 ./ max.(row_measure, eps(Float64))
        A0 = row_scale .* A0
        A1 = row_scale .* A1
        A2 = row_scale .* A2

        column_measure = vec(maximum(abs.(vcat(A0, A1, A2)); dims=1))
        column_scale = 1 ./ max.(column_measure, eps(Float64))
        A0 = A0 .* transpose(column_scale)
        A1 = A1 .* transpose(column_scale)
        A2 = A2 .* transpose(column_scale)
    end

    zero_block = zeros(ComplexF64, m, m)
    identity_block = Matrix{ComplexF64}(I, m, m)
    left = [zero_block identity_block; -A0 -A1]
    right = [identity_block zero_block; zero_block A2]
    decomposition = eigen(left, right)
    candidates = closest_finite(decomposition.values, target)
    results = []
    for index in candidates
        alpha = decomposition.values[index]
        vector = column_scale .* decomposition.vectors[1:m, index]
        push!(results, (alpha=alpha, vector=vector))
    end
    return results
end

function iar_candidates(A0, A1, A2, target; neigs=4, maxit=500, tol=1e-11,
                        initial_vector=nothing)
    problem = PEP([A0, A1, A2])
    initial_vector = isnothing(initial_vector) ? randn(size(A0, 1)) : initial_vector
    values, vectors = iar(
        problem, σ=target, neigs=neigs, maxit=maxit, tol=tol, v=initial_vector,
    )
    order = sortperm(eachindex(values); by=i -> abs(values[i] - target))
    return [(alpha=values[i], vector=vectors[:, i]) for i in order]
end


function run_finite_domain_diagnostics(;
        resolutions=(29, 39, 49, 59, 69, 79, 99, 119, 149, 199),
        target=0.13263552)
    base = prepare_baseflow()
    println("Finite-domain Type-II test in sim coordinate, ymax=20")
    flush(stdout)
    for N in resolutions
        A0, A1, A2, _ = matrices_at_resolution(
            base, N; coordinate="sim", regularization=0.0,
            grid="finite", ymax=20.0,
        )
        matrices = (A0, A1, A2)
        candidates = iar_candidates(
            A0, A1, A2, target; neigs=1, maxit=800, tol=1e-11,
            initial_vector=ones(size(A0, 1)),
        )
        print_candidates("N=$(N), eps=0", candidates, matrices; limit=1)
    end
end

function run_domain_length_diagnostics(;
        resolutions=(59, 99, 149), domains=(20.0, 30.0, 40.0),
        target=0.13263552)
    base = prepare_baseflow()
    println("Finite-domain length test for Type II in sim coordinate")
    flush(stdout)
    for ymax in domains
        println("  ymax=", ymax)
        flush(stdout)
        for N in resolutions
            A0, A1, A2, _ = matrices_at_resolution(
                base, N; coordinate="sim", regularization=0.0,
                grid="finite", ymax=ymax,
            )
            matrices = (A0, A1, A2)
            candidates = iar_candidates(
                A0, A1, A2, target; neigs=1, maxit=800, tol=1e-11,
                initial_vector=ones(size(A0, 1)),
            )
            print_candidates(
                "N=$(N), eps=0", candidates, matrices; limit=1,
            )
        end
    end
end

function run_high_resolution_diagnostics(;
        resolutions=(99, 119, 149, 199), target=0.13263552)
    base = prepare_baseflow()
    println("High-resolution Type-II test in sim coordinate")
    flush(stdout)
    for N in resolutions
        println("  N=", N)
        flush(stdout)
        for epsilon in (0.0, 1e-10, 1e-8, 1e-6)
            A0, A1, A2, clipped = matrices_at_resolution(
                base, N; coordinate="sim", regularization=epsilon,
            )
            matrices = (A0, A1, A2)
            try
                candidates = iar_candidates(
                    A0, A1, A2, target; neigs=1, maxit=800, tol=1e-11,
                )
                print_candidates(
                    "eps=$(epsilon), clipped=$(clipped)", candidates, matrices;
                    limit=1,
                )
            catch error
                println(
                    "    eps=", epsilon, " failed: ", sprint(showerror, error),
                )
                flush(stdout)
            end
        end
    end
end

function print_candidates(label, candidates, matrices; limit=3)
    A0, A1, A2 = matrices
    print("    ", label, ":")
    for result in first(candidates, min(limit, length(candidates)))
        residual = polynomial_residual(A0, A1, A2, result.alpha, result.vector)
        print("  ", result.alpha, " [res=", residual, "]")
    end
    println()
    flush(stdout)
end

function run_diagnostics(; resolutions=(29, 39, 49, 59, 69, 79, 99), target=0.13228)
    base = prepare_baseflow()
    coordinate = "sim"
    println("Coordinate: ", coordinate)
    flush(stdout)
    for N in resolutions
        A0, A1, A2, clipped = matrices_at_resolution(base, N; coordinate=coordinate)
        matrices = (A0, A1, A2)
        println(
            "  N=", N,
            ", dimension=", size(A0, 1),
            ", clipped points=", clipped,
            ", cond(A2)=", cond(A2),
        )
        flush(stdout)
        try
            print_candidates("IAR eps=1e-8", iar_candidates(A0, A1, A2, target), matrices)
        catch error
            println("    IAR failed: ", sprint(showerror, error))
            flush(stdout)
        end

        if N <= 69
            A0s, A1s, A2s, _ = matrices_at_resolution(
                base, N; coordinate=coordinate, regularization=0.0,
            )
            singular_matrices = (A0s, A1s, A2s)
            print_candidates(
                "balanced QZ eps=0",
                dense_companion(A0s, A1s, A2s, target; equilibrate=true),
                singular_matrices,
            )
        end
    end
end

if !isempty(ARGS) && ARGS[1] == "high"
    run_high_resolution_diagnostics()
elseif !isempty(ARGS) && ARGS[1] == "finite"
    run_finite_domain_diagnostics()
elseif !isempty(ARGS) && ARGS[1] == "domain"
    run_domain_length_diagnostics()
else
    run_diagnostics()
end

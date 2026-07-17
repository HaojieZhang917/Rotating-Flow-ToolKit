using LinearAlgebra

include("CRD_STA.jl")

function smooth_test_profiles(N)
    D, D2, y = CRD_BF.Cheb(N)
    yv = vec(y)
    T = reshape(1 .+ 0.15 .* exp.(-yv), :, 1)
    rho = 1 ./ T
    F = reshape(0.2 .* yv .* exp.(-yv), :, 1)
    G = reshape(-1 .+ exp.(-yv), :, 1)
    H = reshape(-0.4 .* (1 .- exp.(-yv)), :, 1)
    sigma = 0.72
    lambda = -(2 / 3) .* T
    kappa = T ./ sigma
    return F, G, H, rho, lambda, kappa, T, D, D2
end

function relative_error(A, B)
    return norm(A - B) / max(norm(A), norm(B), eps(Float64))
end

function spatial_implementation_audit(; N=10)
    F, G, H, rho, lambda, kappa, T, D, D2 = smooth_test_profiles(N)
    sigma = 0.72
    gamma = 1.4
    R = 300.0
    Ma = 0.3 / R
    beta = 0.06
    omega = 0.02
    Ro = -1.0
    Co = 2.0

    direct = Spatial_mode_BEK1(
        F, G, H, rho, lambda, kappa, T, sigma, gamma,
        R, Ma, omega, beta, N, Ro, Co, D, D2,
    )
    coefficients = Spatial_mode_BEK(
        F, G, H, rho, lambda, kappa, T, sigma, gamma,
        R, Ma, N, Ro, Co, D, D2; regularization=0.0,
    )
    assembled = assemble_mat(coefficients, D, D2, beta, omega)

    println("Spatial implementation comparison")
    for (name, A, B) in zip(("L0", "L1", "L2"), direct, assembled)
        delta = A - B
        significant = findall(abs.(delta) .> 1e-12)
        println(
            "  ", name,
            ": relative error = ", relative_error(A, B),
            ", max difference = ", maximum(abs.(delta)),
            ", entries > 1e-12 = ", length(significant),
        )
        if !isempty(significant)
            println("    first differing indices: ", first(significant, min(8, length(significant))))
        end
        n = N + 1
        for row_block in 1:5, column_block in 1:5
            rows = (row_block - 1) * n + 1:row_block * n
            columns = (column_block - 1) * n + 1:column_block * n
            block_error = maximum(abs.(delta[rows, columns]))
            if block_error > 1e-10
                println(
                    "    block (", row_block, ",", column_block,
                    ") max difference = ", block_error,
                )
            end
        end
    end
end

function temporal_implementation_audit(; N=10)
    F, G, H, rho, lambda, kappa, T, D, D2 = smooth_test_profiles(N)
    sigma = 0.72
    gamma = 1.4
    R = 300.0
    Ma = 0.3 / R
    alpha = 0.18
    beta = 0.06
    Ro = -1.0
    Co = 2.0

    actual_B0, actual_B1 = Timemode(
        F, G, H, rho, lambda, kappa, T, sigma, gamma,
        R, Ma, alpha, beta, N, Ro, Co, D, D2,
    )

    coefficients = Spatial_mode_BEK(
        F, G, H, rho, lambda, kappa, T, sigma, gamma,
        R, Ma, N, Ro, Co, D, D2,
    )
    L0, L1, L2 = assemble_mat(coefficients, D, D2, beta, 0.0)
    expected_B0 = L0 + alpha .* L1 + alpha^2 .* L2
    expected_B1 = im .* coefficients.Ta

    unused = zeros(ComplexF64,size(expected_B0))
    expected_B0, expected_B1, _ = boundary_condition(
        expected_B0, expected_B1, unused, N,
    )

    println("Temporal implementation comparison")
    for (name, A, B) in (
        ("B0", actual_B0, expected_B0),
        ("B1", actual_B1, expected_B1),
    )
        delta = A - B
        println(
            "  ", name,
            ": relative error = ", relative_error(A, B),
            ", max difference = ", maximum(abs.(delta)),
            ", entries > 1e-10 = ", count(abs.(delta) .> 1e-10),
        )
    end
end

spatial_implementation_audit()
temporal_implementation_audit()

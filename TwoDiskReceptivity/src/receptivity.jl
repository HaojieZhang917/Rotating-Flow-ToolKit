# Receptivity analysis functions

using LinearAlgebra

"""
    normalise_eigenvectors!(eigvec_direct, eigvec_adjoint, N_cheb)

Normalise the direct and adjoint eigenvectors:

- Direct eigenvector: normalised such that max|v| = 1
- Adjoint eigenvector: normalised such that max|u| = 1

Returns the normalised eigenvectors and their full velocity profiles.
"""
function normalise_eigenvectors!(eigvec_direct, eigvec_adjoint, N_cheb)
    vel_direct = reconstruct_eigenvector(eigvec_direct, N_cheb, 1)
    vel_adjoint = reconstruct_eigenvector(eigvec_adjoint, N_cheb, 1)

    # Normalise direct eigenvector by max of v-component
    scale_direct = findmax(abs.(vel_direct[2]))[1]
    eigvec_direct_norm = eigvec_direct ./ scale_direct

    # Normalise adjoint eigenvector by max of u-component
    scale_adjoint = findmax(abs.(vel_adjoint[1]))[1]
    eigvec_adjoint_norm = eigvec_adjoint ./ scale_adjoint

    # Reconstruct normalised profiles
    vel_direct_norm = reconstruct_eigenvector(eigvec_direct_norm, N_cheb, 1)
    vel_adjoint_norm = reconstruct_eigenvector(eigvec_adjoint_norm, N_cheb, 1)

    return eigvec_direct_norm, eigvec_adjoint_norm, vel_direct_norm, vel_adjoint_norm
end

"""
    chebyshev_quadrature_weights(N::Int)

Compute Clenshaw-Curtis quadrature weights for the Chebyshev grid of size N+1,
then form the block-diagonal mass matrix for the 4-component state vector.

# Returns
- `W::Matrix{Float64}`: Quadrature weight matrix (reduced by boundary conditions)
"""
function chebyshev_quadrature_weights(N::Int)
    N_cheb = N
    x = [cos(π * j / N) for j in 0:N]

    w = zeros(N + 1)
    for j in 0:N
        s = 0.0
        for k in 1:floor(Int, N/2)
            term = 2.0 / (1.0 - (2k)^2) * cos(2 * k * j * π / N)
            if 2k == N
                s += 0.5 * term
            else
                s += term
            end
        end
        c_j = (j == 0 || j == N) ? 1.0 : 2.0
        w[j+1] = (c_j / N) * (1.0 + s)
    end

    # Jacobian for mapping to physical domain (mode 1: [0,1] → dx/dξ = 1/2)
    Jac = zeros(N_cheb + 1, N_cheb + 1)
    x = reverse(x)
    for i in 1:N_cheb+1
        Jac[i, i] = 1/2
    end

    W = kron(I(4), diagm(Jac * w))

    # Remove boundary rows/columns
    bc_rows = (1, N_cheb + 1, N_cheb + 2, 2N_cheb + 2, 2N_cheb + 3, 3N_cheb + 3, 4N_cheb + 4)
    interior = setdiff(1:4*(N_cheb+1), bc_rows)
    W = W[interior, interior]

    return W
end

"""
    compute_coupling_coefficient(eigvec_adjoint, eigvec_direct, eigval_direct, eigval_adjoint, L1, L2, W)

Compute the receptivity coupling coefficient Q:

    Q = φ† · W · (L1 + (α + α†) L2) · φ

where φ is the direct eigenvector and φ† the adjoint eigenvector.
"""
function compute_coupling_coefficient(eigvec_adjoint, eigvec_direct, eigval_direct, eigval_adjoint, L1, L2, W)
    M = L1 + (eigval_direct[1] + eigval_adjoint[1]) * L2
    Q = transpose(eigvec_adjoint[:, 1]) * W * M * eigvec_direct[:, 1]
    return Q[1]
end

"""
    compute_receptivity_coefficient(params, hr, ls, eigval_direct, F, G, R, vel_adjoint, Q)

Compute the receptivity coefficient Cr from wall roughness.

The receptivity coefficient quantifies how wall roughness (height hr, length scale ls)
couples to the instability mode.

# Arguments
- `params::FlowParameters`: Flow parameters
- `hr`: Roughness height
- `ls`: Roughness length scale (Gaussian width)
- `eigval_direct`: Direct eigenvalue array
- `F, G`: Base flow velocity profiles
- `R`: Reynolds number
- `vel_adjoint`: Normalised adjoint velocity profile (u†, v†, w†, p†)
- `Q`: Coupling coefficient

# Returns
- `Cr`: Receptivity coefficient (real, positive)
- `BC`: Boundary condition contribution (complex)
"""
function compute_receptivity_coefficient(params::FlowParameters, hr, ls, eigval_direct, F, G, D, Res, vel_adjoint, Q)
    # Gaussian roughness spectrum
    Hx = hr * exp(-(eigval_direct[1])^2 / (4*ls)) * sqrt(π / ls)

    # Wall velocity perturbations from roughness
    u_wall = -(D * F)[1] * Hx
    v_wall = -(D * G)[1] * Hx

    # Adjoint boundary terms
    item_2 = (1 / params.R) * (D * vel_adjoint[1])[1] * u_wall
    item_3 = (1 / params.R) * (D * vel_adjoint[2])[1] * v_wall
    BC = (item_2 + item_3) / sqrt(Res)

    Cr = abs(-im * BC / Q)

    return Cr, BC
end

"""
    solve_receptivity(params::FlowParameters, baseflow, grid, sigma, hr, ls)

Complete receptivity analysis pipeline.

# Returns
- `ReceptivityResult`
"""
function solve_receptivity(params::FlowParameters, baseflow, grid, sigma, hr, ls)
    D = grid.D

    # Step 1: Solve direct + adjoint eigenvalue problems
    eigval_dir, eigvec_dir, eigval_adj, eigvec_adj, L1, L2 =
        solve_stability_eigenproblem(params, baseflow, grid, sigma)

    # Step 2: Normalise eigenvectors
    eigvec_dir_norm, eigvec_adj_norm, vel_dir, vel_adj =
        normalise_eigenvectors!(eigvec_dir, eigvec_adj, params.N_cheb)

    # Step 3: Compute quadrature weights
    W = chebyshev_quadrature_weights(params.N_cheb)

    # Step 4: Coupling coefficient Q
    Q = compute_coupling_coefficient(eigvec_adj_norm, eigvec_dir_norm,
                                      eigval_dir, eigval_adj, L1, L2, W)

    # Step 5: Receptivity coefficient Cr
    Cr, BC = compute_receptivity_coefficient(params, hr, ls, eigval_dir,
                                              baseflow.F, baseflow.G, D,
                                              params.Res, vel_adj, Q)

    return ReceptivityResult(
        eigval_dir[1], eigval_adj[1],
        eigvec_dir[:, 1], eigvec_adj[:, 1],
        Cr, BC, Q,
        vel_dir, vel_adj
    )
end

# Eigenvalue problem assembly and solution

using LinearAlgebra
using NonlinearEigenproblems

"""
    apply_boundary_conditions!(L0, L1, L2, N_cheb)

Apply homogeneous Dirichlet boundary conditions by removing rows and columns
corresponding to boundary nodes from the PEP matrices.

The boundary conditions are: u=v=w=0 at both walls (z=0 and z=1).
"""
function apply_boundary_conditions!(L0, L1, L2, N_cheb)
    bc_rows = (1, N_cheb + 1, N_cheb + 2, 2N_cheb + 2, 2N_cheb + 3, 3N_cheb + 3, 4N_cheb + 4)
    interior = setdiff(1:4*(N_cheb+1), bc_rows)

    L0 = L0[interior, interior]
    L1 = L1[interior, interior]
    L2 = L2[interior, interior]

    return L0, L1, L2
end

"""
    reconstruct_eigenvector(eigvec_reduced, N_cheb, k)

Reconstruct the full 4-component eigenvector (u, v, w, p) on the physical grid
from the reduced eigenvector (after boundary condition removal).

# Arguments
- `eigvec_reduced`: Eigenvector in reduced space (size 4N-6)
- `N_cheb`: Number of Chebyshev points
- `k`: Which eigenpair to extract

# Returns
- `(u, v, w, p)`: Four velocity/pressure components on the full grid (length N+1 each)
"""
function reconstruct_eigenvector(eigvec_reduced, N_cheb, k)
    N = N_cheb + 1
    v = eigvec_reduced[:, k]

    u = zeros(ComplexF64, N)
    v_vel = zeros(ComplexF64, N)
    w = zeros(ComplexF64, N)
    p = zeros(ComplexF64, N)

    # Interior points: the reduced vector is laid out as [u_interior; v_interior; w_interior; p_interior]
    u[2:N-1] = v[1:N-2]
    v_vel[2:N-1] = v[N-1:2N-4]
    w[2:N-1] = v[2N-3:3N-6]
    p[1:N-1] = v[3N-5:4N-7]

    return (u, v_vel, w, p)
end

"""
    solve_eigenvalue_problem(L0, L1, L2, sigma, nev)

Solve the quadratic eigenvalue problem using the Infinite Arnoldi Method (IAR)
with a shift-and-invert strategy.

# Arguments
- `L0, L1, L2`: PEP coefficient matrices (after BCs applied)
- `sigma`: Shift-invert target (complex)
- `nev`: Number of eigenvalues to compute

# Returns
- `eigvals, eigvecs`: Eigenvalues and eigenvectors
"""
function solve_eigenvalue_problem(L0, L1, L2, sigma::ComplexF64, nev::Int)
    nep = PEP([L0, L1, L2])
    eigvals, eigvecs = iar(nep; σ=sigma, neigs=nev, maxit=500, tol=1e-12)
    return eigvals, eigvecs
end

"""
    solve_stability_eigenproblem(params::FlowParameters, baseflow::BaseFlowProfile, grid::ChebyshevGrid, sigma)

Solve both the direct and adjoint eigenvalue problems for the stability analysis.

# Returns
- `eigval_direct`, `eigvec_direct`: Direct eigenvalues/vectors
- `eigval_adjoint`, `eigvec_adjoint`: Adjoint eigenvalues/vectors
- `L1, L2`: PEP matrices L1, L2 (used in receptivity coupling coefficient)
"""
function solve_stability_eigenproblem(params::FlowParameters, baseflow::BaseFlowProfile, grid::ChebyshevGrid, sigma::ComplexF64)
    F, G, H = baseflow.F, baseflow.G, baseflow.H
    D, D2 = grid.D, grid.D2
    R = params.R
    N_cheb = params.N_cheb
    Res = params.Res
    Ro = params.Ro
    omega = params.OMEGA / R
    be = params.n / R

    # Assemble coefficient matrices
    Co = 2 - Ro - Ro^2
    if params.mode == 1
        cof = assemble_coeff_matrix_BEK1(F, G, H, R, N_cheb, D, D2, Res)
    else
        cof = assemble_coeff_matrix_BEK2(F, G, H, R, N_cheb, D, D2, Ro, Co)
    end

    # Direct problem
    L0, L1, L2 = assemble_direct_matrices(cof, D, D2, be, omega, R)
    L0, L1, L2 = apply_boundary_conditions!(L0, L1, L2, N_cheb)
    eigval_direct, eigvec_direct = solve_eigenvalue_problem(L0, L1, L2, sigma, 1)

    # Adjoint problem
    A0, A1, A2 = assemble_adjoint_matrices(cof, D, D2, be, omega, R)
    A0, A1, A2 = apply_boundary_conditions!(A0, A1, A2, N_cheb)
    eigval_adjoint, eigvec_adjoint = solve_eigenvalue_problem(A0, A1, A2, sigma, 1)

    return eigval_direct, eigvec_direct, eigval_adjoint, eigvec_adjoint, L1, L2
end

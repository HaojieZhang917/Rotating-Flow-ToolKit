using TwoDiskReceptivity
using Test
using LinearAlgebra

@testset "TwoDiskReceptivity.jl" begin

    @testset "Base flow computation" begin
        params = FlowParameters(1000, -1.0, 0.0, 300, 30, 0.0, 129, 1)
        baseflow, grid = solve_baseflow(params)

        @test size(baseflow.F) == (130, 1)
        @test size(baseflow.G) == (130, 1)
        @test size(baseflow.H) == (130, 1)
        @test length(grid.z) == 130
        @test size(grid.D) == (130, 130)
        @test size(grid.D2) == (130, 130)

        # Boundary values: G(z=0) ≈ 0 after G.-1 adjustment (Bödewadt condition)
        @test abs(baseflow.G[1, 1]) < 0.01
    end

    @testset "Coefficient matrix assembly" begin
        params = FlowParameters(1000, -1.0, 0.0, 300, 30, 0.0, 129, 1)
        baseflow, grid = solve_baseflow(params)
        F, G, H = baseflow.F, baseflow.G, baseflow.H
        D, D2 = grid.D, grid.D2

        cof = assemble_coeff_matrix_BEK1(F, G, H, params.R, params.N_cheb, D, D2, params.Res)
        full_size = 4 * (params.N_cheb + 1)
        @test size(cof.Ta) == (full_size, full_size)
        @test size(cof.A)  == (full_size, full_size)
        @test size(cof.Vxx) == (full_size, full_size)
    end

    @testset "PEP matrix assembly and boundary conditions" begin
        params = FlowParameters(1000, -1.0, 0.0, 300, 30, 0.0, 129, 1)
        baseflow, grid = solve_baseflow(params)
        F, G, H = baseflow.F, baseflow.G, baseflow.H
        D, D2 = grid.D, grid.D2
        omega = params.OMEGA / params.R
        be = params.n / params.R

        cof = assemble_coeff_matrix_BEK1(F, G, H, params.R, params.N_cheb, D, D2, params.Res)
        L0, L1, L2 = assemble_direct_matrices(cof, D, D2, be, omega, params.R)
        A0, A1, A2 = assemble_adjoint_matrices(cof, D, D2, be, omega, params.R)

        # Before BC: full 4(N+1) × 4(N+1)
        full_size = 4 * (params.N_cheb + 1)
        @test size(L0) == (full_size, full_size)
        @test size(A0) == (full_size, full_size)

        L0, L1, L2 = apply_boundary_conditions!(L0, L1, L2, params.N_cheb)
        A0, A1, A2 = apply_boundary_conditions!(A0, A1, A2, params.N_cheb)

        # After BC: 4(N+1) - 7 (u=v=w=p=0 at both walls)
        expected_size = full_size - 7
        @test size(L0) == (expected_size, expected_size)
        @test size(L1) == (expected_size, expected_size)
        @test size(L2) == (expected_size, expected_size)
        @test size(A0) == (expected_size, expected_size)
        @test size(A1) == (expected_size, expected_size)
        @test size(A2) == (expected_size, expected_size)
    end

    @testset "Adjoint–direct eigenvalue equivalence" begin
        # The adjoint PEP A(α) = L(α)^T should have exactly the same
        # eigenvalues as the direct PEP L(α), since det(A) = det(L^T) = det(L).
        #
        # We verify this by solving both problems independently with IAR
        # and comparing the leading eigenvalues.

        params = FlowParameters(1000, -1.0, 0.0, 300, 30, 0.0, 129, 1)
        baseflow, grid = solve_baseflow(params)
        F, G, H = baseflow.F, baseflow.G, baseflow.H
        D, D2 = grid.D, grid.D2
        omega = params.OMEGA / params.R
        be = params.n / params.R

        # Assemble coefficient matrices
        cof = assemble_coeff_matrix_BEK1(F, G, H, params.R, params.N_cheb, D, D2, params.Res)

        # Direct PEP
        L0, L1, L2 = assemble_direct_matrices(cof, D, D2, be, omega, params.R)
        L0, L1, L2 = apply_boundary_conditions!(L0, L1, L2, params.N_cheb)

        # Adjoint PEP
        A0, A1, A2 = assemble_adjoint_matrices(cof, D, D2, be, omega, params.R)
        A0, A1, A2 = apply_boundary_conditions!(A0, A1, A2, params.N_cheb)

        sigma = 0.4 + 0.0im
        nev = 2

        eigval_dir, _ = solve_eigenvalue_problem(L0, L1, L2, sigma, nev)
        eigval_adj, _ = solve_eigenvalue_problem(A0, A1, A2, sigma, nev)

        @test length(eigval_dir) >= 1
        @test length(eigval_adj) >= 1

        # Leading eigenvalue should match (within IAR tolerance)
        @test eigval_dir[1] ≈ eigval_adj[1] atol=1e-10 rtol=1e-6

        # If we got two eigenvalues, the second should also match
        if length(eigval_dir) >= 2 && length(eigval_adj) >= 2
            @test eigval_dir[2] ≈ eigval_adj[2] atol=1e-10 rtol=1e-6
        end
    end

    @testset "Eigenvalue problem solution (high-level API)" begin
        params = FlowParameters(1000, -1.0, 0.0, 300, 30, 0.0, 129, 1)
        baseflow, grid = solve_baseflow(params)
        sigma = 0.4 + 0.0im

        eigval_dir, eigvec_dir, eigval_adj, eigvec_adj, _, _ =
            solve_stability_eigenproblem(params, baseflow, grid, sigma)

        @test length(eigval_dir) == 1
        @test size(eigvec_dir, 2) == 1
        @test length(eigval_adj) == 1

        # Direct and adjoint eigenvalues should match
        @test eigval_dir[1] ≈ eigval_adj[1] atol=1e-10 rtol=1e-6
    end

    @testset "Receptivity coefficient" begin
        params = FlowParameters(1000, -1.0, 0.0, 300, 30, 0.0, 129, 1)
        baseflow, grid = solve_baseflow(params)
        sigma = 0.4 + 0.0im
        hr = 1.0 / params.Res
        ls = 0.5

        result = solve_receptivity(params, baseflow, grid, sigma, hr, ls)

        @test result.Cr >= 0.0
        @test length(result.velocity_direct) == 4
        @test length(result.velocity_adjoint) == 4

        # Direct and adjoint eigenvalues in the result should also match
        @test result.eigval_direct ≈ result.eigval_adjoint atol=1e-10 rtol=1e-6
    end

    @testset "Quadrature weights" begin
        W = chebyshev_quadrature_weights(129)
        expected_size = 4*130 - 7
        @test size(W) == (expected_size, expected_size)
        # Weights should be positive
        @test all(diag(W) .>= 0)
    end

    @testset "Eigenvector reconstruction" begin
        # Create a dummy reduced eigenvector
        N_cheb = 129
        reduced_size = 4*(N_cheb+1) - 7
        dummy_eigvec = randn(ComplexF64, reduced_size, 2)

        u, v, w, p = reconstruct_eigenvector(dummy_eigvec, N_cheb, 1)
        @test length(u) == N_cheb + 1
        @test length(v) == N_cheb + 1
        @test length(w) == N_cheb + 1
        @test length(p) == N_cheb + 1
        # Boundary values should be zero
        @test u[1] == 0.0 && u[end] == 0.0
        @test v[1] == 0.0 && v[end] == 0.0
        @test w[1] == 0.0 && w[end] == 0.0
    end

end

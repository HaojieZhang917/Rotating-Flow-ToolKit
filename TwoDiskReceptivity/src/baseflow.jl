# Base flow computation using SciPy's BVP solver

using PyCall
using LinearAlgebra
using BSplineKit

"""
    compute_baseflow(Res::Int, Ro::Float64, Ts::Float64, mode::Int)

Solve the Bödewadt–von Kármán similarity ODEs for the base flow between
two rotating disks using SciPy's `solve_bvp`.

# Arguments
- `Res`: Reynolds number based on disk spacing (√Res appears in scaling)
- `Ro`: Rotation ratio between the two disks
- `Ts`: Suction/injection parameter (wall-normal velocity at the disk)
- `mode`: Base flow mode (1 = single-disk Bödewadt, 2 = two-disk, 3 = two-disk with kappa)

# Returns
- `u0, v0, w0, du0, dv0, x`: Velocity profiles and spatial coordinate from the BVP solution
"""
function compute_baseflow(Res::Int, Ro::Float64, Ts::Float64, mode::Int)
    if mode == 1
        py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        from math import sqrt
        Re_s = $Res
        Ts = $Ts
        def oneDiskODE(z, y):
                # Y0 = H, Y1 = F, Y2 = F', Y3 = F'', Y4 = G, Y5 = G'
                dH = -2 * sqrt(Re_s) * y[1]
                dydz = np.zeros((6, len(z)))
                dydz = np.array([dH, y[2], y[3],
                    Re_s * ((1/sqrt(Re_s)) * (y[3] * y[0] + y[2] * dH) - 2 * y[4] * y[5] + 2 * y[1] * y[2]),
                    y[5],
                    Re_s * ((1/sqrt(Re_s)) * y[5] * y[0] + 2 * y[1] * y[4])])
                return dydz

        def oneDiskBC(ya, yb):
                resa = np.array([ya[0] + Ts, ya[1], ya[4] - 1])
                resb = np.array([yb[0], yb[1], yb[4] - 0])
                return np.concatenate((resa, resb))

        z = np.linspace(0, 1, 1000)
        y_guess = np.zeros((6, z.size))
        y_guess[0] = 1
        y_guess[4] = 1
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess, max_nodes=5000000)
        x_plot = np.linspace(0, 1, 2000)
        y1_plot = solution.sol(x_plot)[0]
        y2_plot = solution.sol(x_plot)[1]
        y3_plot = solution.sol(x_plot)[4]
        y4_plot = solution.sol(x_plot)[2]
        y5_plot = solution.sol(x_plot)[5]
        """
        w0  = py"y1_plot"
        u0  = py"y2_plot"
        v0  = py"y3_plot"
        du0 = py"y4_plot"
        dv0 = py"y5_plot"
        x   = py"x_plot"

    elseif mode == 2
        py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        kappa = $Ro
        Ts = $Ts
        def oneDiskODE(z, y):
                # Y0 = H, Y1 = F', Y2 = F, Y3 = G', Y4 = G
                dydz = np.zeros((5, len(z)))
                dydz = np.array([-2*y[2], y[2]*y[2] - y[4]*y[4] + y[0]*y[1],
                                 y[1], 2*y[2]*y[4] + y[3]*y[0], y[3]])
                return dydz

        def oneDiskBC(ya, yb):
                resa = np.array([ya[0]+Ts, ya[2], ya[4]-1.0])
                resb = np.array([yb[2], yb[4]])
                return np.concatenate((resa, resb))

        z = np.linspace(0, 30, 2000)
        y_guess = np.zeros((5, z.size))
        y_guess[0] = 1.2
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess, tol=1e-10, max_nodes=5000000)
        x_plot = np.linspace(0, 30, 2000)
        y1_plot = solution.sol(x_plot)[0]
        y2_plot = solution.sol(x_plot)[2]
        y3_plot = solution.sol(x_plot)[4]
        y4_plot = solution.sol(x_plot)[1]
        y5_plot = solution.sol(x_plot)[3]
        """
        w0  = py"y1_plot"
        u0  = py"y2_plot"
        v0  = py"y3_plot"
        du0 = py"y4_plot"
        dv0 = py"y5_plot"
        x   = py"x_plot"

    elseif mode == 3
        py"""
        import numpy as np
        from scipy.integrate import solve_bvp
        kappa = $Ro
        Ts = $Ts
        def oneDiskODE(z, y):
                # Y0 = H, Y1 = F', Y2 = F, Y3 = G', Y4 = G
                dydz = np.zeros((5, len(z)))
                dydz = np.array([-2.0*y[2],
                    kappa*(y[2]*y[2] + y[0]*y[1] - (y[4]*y[4] - 1.0)) - (2.0 - kappa - kappa**2)*(y[4] - 1.0),
                    y[1],
                    kappa*(2.0*y[2]*y[4] + y[0]*y[3]) + (2.0 - kappa - kappa**2)*y[2],
                    y[3]])
                return dydz

        def oneDiskBC(ya, yb):
                resa = np.array([ya[0]+Ts, ya[2], ya[4]])
                resb = np.array([yb[2], yb[4] - 1.0])
                return np.concatenate((resa, resb))

        z = np.linspace(0, 30, 2000)
        y_guess = np.zeros((5, z.size))
        y_guess[0] = 1.2
        y_guess[4] = 1
        solution = solve_bvp(oneDiskODE, oneDiskBC, z, y_guess, tol=1e-10, max_nodes=5000000)
        x_plot = np.linspace(0, 30, 2000)
        y1_plot = solution.sol(x_plot)[0]
        y2_plot = solution.sol(x_plot)[2]
        y3_plot = solution.sol(x_plot)[4]
        y4_plot = solution.sol(x_plot)[1]
        y5_plot = solution.sol(x_plot)[3]
        """
        w0  = py"y1_plot"
        u0  = py"y2_plot"
        v0  = py"y3_plot"
        du0 = py"y4_plot"
        dv0 = py"y5_plot"
        x   = py"x_plot"
    end
    return u0, v0, w0, du0, dv0, x
end

"""
    chebyshev_grid(N::Int, mode::Int)

Generate a Chebyshev collocation grid and differentiation matrices.

For mode=1, the grid is mapped to [0, 1] (single-disk Bödewadt).
For other modes, a stretching transformation maps the grid to the semi-infinite domain.

# Returns
- `ChebyshevGrid` containing z, D, D2
"""
function chebyshev_grid(N::Int, mode::Int)
    θ = range(0, length=N+1, stop=π)
    x = reshape(-cos.(θ), N+1, 1)
    c = [2; ones(N-1, 1); 2] .* (-1) .^ (0:N)
    X = repeat(x, 1, N+1)
    dX = X - X'
    D = (c * (1 ./ c)') ./ (dX .+ I(N+1))
    D = D - diagm(vec(sum(D, dims=2)))

    if mode == 1
        D = 2 * D
        x = 0.5 * (x .+ 1)
        D2 = D^2
    else
        a = 2.0
        b = 0.6
        c0 = 0.5
        for i = 1:N+1
            D[i, :] = D[i, :] .* (1 - b*x[i] - (1-b)*(x[i]^3 + c0*(1-x[i]^2)))^2 /
                      (2a * (b + 3*(1-b)*x[i]^2 - 2*c0*(1-b)*x[i]))
        end
        for i = 1:N+1
            x[i] = a * (1 + b*x[i] + (1-b)*(x[i]^3 + c0*(1-x[i]^2))) /
                   (1 - b*x[i] - (1-b)*(x[i]^3 + c0*(1-x[i]^2)))
            if x[i] > 30
                x[i] = 30
            end
        end
        D2 = D^2
    end
    return ChebyshevGrid(vec(x), D, D2)
end

"""
    interpolate_to_cheb(u0, v0, w0, x_bvp, z_cheb::Vector{Float64}, N::Int)

Interpolate the BVP solution onto the Chebyshev grid using B-splines.

# Returns
- `BaseFlowProfile` containing F, G, H on the Chebyshev grid
"""
function interpolate_to_cheb(u0, v0, w0, x_bvp, z_cheb::Vector{Float64}, N::Int)
    F = zeros(N+1, 1)
    G = zeros(N+1, 1)
    H = zeros(N+1, 1)

    x = range(minimum(x_bvp), maximum(x_bvp), length=length(x_bvp))

    it_u = BSplineKit.interpolate(x, u0, BSplineOrder(4))
    it_v = BSplineKit.interpolate(x, v0, BSplineOrder(4))
    it_w = BSplineKit.interpolate(x, w0, BSplineOrder(4))

    for i = 1:N+1
        F[i, 1] = it_u(z_cheb[i])
        G[i, 1] = it_v(z_cheb[i])
        H[i, 1] = it_w(z_cheb[i])
    end

    return BaseFlowProfile(F, G .- 1, H)
end

"""
    solve_baseflow(params::FlowParameters)

Convenience function that computes the base flow and returns both the profile
and the Chebyshev grid.

# Returns
- `baseflow::BaseFlowProfile`, `grid::ChebyshevGrid`
"""
function solve_baseflow(params::FlowParameters)
    u0, v0, w0, _, _, x_bvp = compute_baseflow(params.Res, params.Ro, params.Ts, params.mode)
    grid = chebyshev_grid(params.N_cheb, params.mode)
    baseflow = interpolate_to_cheb(u0, v0, w0, x_bvp, grid.z, params.N_cheb)
    return baseflow, grid
end

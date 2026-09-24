module CompressRotatingCavity

using LinearAlgebra
using BSplineKit
using NonlinearEigenproblems

include("CRC_BF.jl")
include("CRC_STA.jl")
include("TraditionalBoussinesq.jl")
include("BlackburnBoussinesq.jl")
include("TraditionalBoussinesqTemporal.jl")
include("TraditionalBoussinesqLST.jl")

using .TraditionalBoussinesq: ThermalConfig, ThermalProfile,
    isothermal_profile, fixed_pressure_profile, pressure_tangent,
    fixed_temperature_profile, trace_principal_fold
using .TraditionalBoussinesqTemporal: temporal_operator,
    leading_temporal_eigenvalue
using .BlackburnBoussinesq: BlackburnConfig, BlackburnProfile
using .TraditionalBoussinesqLST: lst_coefficients, spatial_matrices,
    temporal_matrices, expand_mode

const blackburn_isothermal_profile = BlackburnBoussinesq.isothermal_profile
const blackburn_fixed_pressure_profile = BlackburnBoussinesq.fixed_pressure_profile
const blackburn_fixed_temperature_profile = BlackburnBoussinesq.fixed_temperature_profile
const blackburn_pressure_tangent = BlackburnBoussinesq.pressure_tangent
const trace_blackburn_fold = BlackburnBoussinesq.trace_principal_fold

export baseflow, eigsolve, normalize_modes, cheb_quad,
       q_receptivity, cr_receptivity, benchmark_case,
       ThermalConfig, ThermalProfile, isothermal_profile,
       fixed_pressure_profile, fixed_temperature_profile, pressure_tangent,
       trace_principal_fold,
       BlackburnConfig, BlackburnProfile,
       blackburn_isothermal_profile, blackburn_fixed_pressure_profile,
       blackburn_fixed_temperature_profile,
       blackburn_pressure_tangent, trace_blackburn_fold,
       temporal_operator, leading_temporal_eigenvalue,
       lst_coefficients, spatial_matrices, temporal_matrices, expand_mode

"""Compute the legacy isothermal base flow on the Chebyshev grid."""
function baseflow(Res, Ro, Ts, mode, N)
    u0, v0, w0, du0, dv0, x = CRC_BF.BaseFlow(Res, Ro, Ts, mode)
    D, D2, z = CRC_BF.Cheb(N, mode)
    F, G, H = CRC_BF.interp(u0, v0, w0, z, N, mode)
    return (; F, G, H, D, D2, z)
end

"""Solve the quadratic spatial eigenproblem and its transpose adjoint."""
function eigsolve(flow, Res, Ro, mode, c, R, omega, beta;
                  maxit=500, tol=1e-14)
    (; F, G, H, D, D2, z) = flow
    N = length(z) - 1
    cof = CRC_STA.Spatial_mode_BEK1(F, G .- 1, H, R, N, D, D2, Res)

    L0_raw, L1_raw, L2_raw = CRC_STA.assemble_mat(cof, D, D2, beta, omega, R)
    L0, L1, L2 = CRC_STA.boudary_condition(L0_raw, L1_raw, L2_raw, N, mode)
    eigval, eigvec = iar(PEP([L0, L1, L2]); σ=c, neigs=1,
                          maxit=maxit, tol=tol)

    # This is the adjoint assembly used inside notebook cell 1.  The notebook
    # cell 15 helper has a different historical scaling for two terms, so the
    # benchmark keeps the cell-1 form explicitly.
    A0_raw = transpose(cof.D1) + (im * R * beta * transpose(cof.B)) -
        (im * omega * transpose(cof.Ta)) -
        (beta^2 * R^2 * transpose(cof.Vyy)) - transpose(cof.dC) -
        (im * beta * transpose(cof.dVyz)) + transpose(cof.d2Vzz) -
        (transpose(cof.C) + im * beta * R * transpose(cof.Vyz) -
         2 * transpose(cof.dVzz)) * kron(I(4), D) +
        transpose(cof.Vzz) * kron(I(4), D2)
    A1_raw = (im * transpose(cof.A)) - (beta * R * transpose(cof.Vxy)) -
        (im * transpose(cof.dVxz)) - (im * transpose(cof.Vxz)) * kron(I(4), D)
    A2_raw = -transpose(cof.Vxx)
    A0, A1, A2 = CRC_STA.boudary_condition(A0_raw, A1_raw, A2_raw, N, mode)
    eigval_A, eigvec_A = iar(PEP([A0, A1, A2]); σ=c, neigs=1,
                              maxit=maxit, tol=tol)
    return (; eigval, eigval_A, eigvec, eigvec_A, L=(L1, L2), z)
end

function normalize_modes(solution)
    (; eigvec, eigvec_A, z) = solution
    N = length(z) - 1
    vel = CRC_STA.eig_full(eigvec, N, 1)
    vel_A = CRC_STA.eig_full(eigvec_A, N, 1)
    eigvec = eigvec ./ maximum(abs, vel[2])
    eigvec_A = eigvec_A ./ maximum(abs, vel_A[1])
    vel = CRC_STA.eig_full(eigvec, N, 1)
    vel_A = CRC_STA.eig_full(eigvec_A, N, 1)
    return (; solution..., eigvec, eigvec_A, vel, vel_A)
end

"""Chebyshev quadrature matrix with the boundary degrees removed."""
function cheb_quad(N::Int)
    w = zeros(N + 1)
    for j in 0:N
        s = 0.0
        for k in 1:floor(Int, N / 2)
            term = 2.0 / (1.0 - (2k)^2) * cos(2k * j * pi / N)
            s += (2k == N) ? 0.5 * term : term
        end
        c_j = (j == 0 || j == N) ? 1.0 : 2.0
        w[j + 1] = (c_j / N) * (1.0 + s)
    end
    W = kron(I(4), Diagonal(0.5 .* w))
    removed = (1, N + 1, N + 2, 2N + 2, 2N + 3, 3N + 3, 4N + 4)
    keep = setdiff(1:size(W, 1), removed)
    return W[keep, keep]
end

function q_receptivity(solution, W)
    (; eigvec, eigvec_A, eigval, eigval_A, L) = solution
    return transpose(eigvec_A[:, 1]) * W *
           (L[1] + (eigval[1] + eigval_A[1]) * L[2]) * eigvec[:, 1]
end

function cr_receptivity(solution, flow, Res, R, hr, ls, Q)
    (; eigval, vel_A) = solution
    (; F, G, D) = flow
    Hx = hr * exp(-(eigval[1])^2 / (4ls)) * sqrt(pi / ls)
    u_wall = -(D * F)[1] * Hx
    v_wall = -(D * G)[1] * Hx
    BC = ((D * vel_A[1])[1] * u_wall +
          (D * vel_A[2])[1] * v_wall) / (R * sqrt(Res))
    return (; Cr=abs(-im * BC / Q[1]), BC)
end

function benchmark_case(; Res=2000, N=129, Ro=-1.0, R=330.0,
                         beta=0.0503, OMEGA=8.0, Ts=0.0, c=0.5,
                         mode=1, hr=1 / Res, ls=0.5)
    flow = (; baseflow(Res, Ro, Ts, mode, N)..., Res)
    solution = normalize_modes(eigsolve(flow, Res, Ro, mode, c, R,
                                        OMEGA / R, beta))
    Q = q_receptivity(solution, cheb_quad(N))
    receptivity = cr_receptivity(solution, flow, Res, R, hr, ls, Q)
    return (; flow, solution, Q, receptivity)
end

end

module CRC_RA

using LinearAlgebra

export boundary_indices, interior_indices, wall_boundary_vector,
       lifted_response, spatial_residue

"""Indices removed by the legacy four-field local operator.

The ordering is (u,v,w,p), each field having N+1 Chebyshev points.  The
velocity boundary values are prescribed at z=0 and z=1; the last pressure
degree of freedom fixes the pressure gauge.
"""
function boundary_indices(N::Integer)
    n = N + 1
    return [1, n, n + 1, 2n, 2n + 1, 3n, 4n]
end

function interior_indices(N::Integer)
    all = collect(1:4(N + 1))
    b = boundary_indices(N)
    return setdiff(all, b)
end

"""Build boundary data in the ordering used by `boundary_indices`.

`uw`, `vw`, `ww` are complex Fourier amplitudes at the rotor wall.  The
stator values are zero and the final pressure entry is the gauge value.
"""
function wall_boundary_vector(N::Integer; uw=0.0 + 0.0im,
                             vw=0.0 + 0.0im,
                             ww=0.0 + 0.0im,
                             gauge=0.0 + 0.0im)
    return ComplexF64[uw, 0.0 + 0.0im, vw, 0.0 + 0.0im,
                     ww, 0.0 + 0.0im, gauge]
end

"""Solve the inhomogeneous boundary problem by a lifting/partition.

For a raw collocation matrix P whose boundary rows are not used by the
homogeneous stability problem, impose q_B=b and solve the retained interior
equations P_II q_I = -P_IB b.  The result is independent of the particular
interior shape chosen for a lifting function.
"""
function lifted_response(P::AbstractMatrix, N::Integer, b::AbstractVector)
    bi = boundary_indices(N)
    length(b) == length(bi) || throw(DimensionMismatch("boundary vector"))
    ii = interior_indices(N)
    q = zeros(ComplexF64, size(P, 1))
    q[bi] .= ComplexF64.(b)
    Pii = P[ii, ii]
    Pib = P[ii, bi]
    q[ii] .= -(Pii \ (Pib * q[bi]))
    return q
end

"""Extract the discrete spatial-pole receptivity coefficient.

The current paper uses the bilinear transpose convention for the spatial
adjoint.  `v` is the reduced direct eigenvector associated with `alpha` and
`b` is the physical wall-forcing vector at the pole.  The coefficient is
`-im*Bw/Q` for the Fourier/residue convention used in the manuscript.
"""
function spatial_residue(L0::AbstractMatrix, L1::AbstractMatrix,
                         L2::AbstractMatrix, N::Integer,
                         alpha, v::AbstractVector, b::AbstractVector)
    ii = interior_indices(N)
    bi = boundary_indices(N)
    length(v) == length(ii) || throw(DimensionMismatch("direct vector"))
    length(b) == length(bi) || throw(DimensionMismatch("boundary vector"))
    P = L0 + alpha .* L1 + alpha^2 .* L2
    Pprime = L1 + 2alpha .* L2
    Pii = P[ii, ii]
    Pib = P[ii, bi]
    Ppii = Pprime[ii, ii]
    # Algebraic left null vector: transpose is intentional and matches the
    # spatial bilinear adjoint convention used by the existing benchmark.
    left = nullspace(transpose(Pii))
    size(left, 2) == 1 || throw(ArgumentError("eigenvalue is not simple"))
    ell = left[:, 1]
    Q = transpose(ell) * Ppii * v
    Bw = transpose(ell) * Pib * ComplexF64.(b)
    return (; Cr=-im * Bw[1] / Q[1], Bw=Bw[1], Q=Q[1],
            response=lifted_response(P, N, b), left=ell)
end

end

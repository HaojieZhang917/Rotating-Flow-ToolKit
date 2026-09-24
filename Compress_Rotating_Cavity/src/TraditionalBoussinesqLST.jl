module TraditionalBoussinesqLST

using LinearAlgebra

export LSTCoefficients, lst_coefficients, spatial_matrices, temporal_matrices,
       expand_mode

"""
Local parallel-flow LST for the *traditional* centrifugal Boussinesq closure.

The finite-gap base flow is `(H,F,G,T)` on `z/h ∈ [0,1]`, with
`Re_h = Ωh²/ν`, `δ = h/√Re_h`, and `R = r₀/δ`. The local radial coordinate is
`x = (r-r₀)/δ`; a mode is `q(z)exp(iαx+iRβφ-iωt*)`, where `t* = ΩRt`.
Velocity disturbances use `Ωδ`, temperature uses the reference temperature,
and pressure uses `ρΩ²δ²R`. The base flow is already solved and contains no
disturbance wavenumber. Only this LST module receives `α`, `β`, or `ω`.

The five fields are `(u,v,w,θ,p)`. At both walls `u=v=w=θ=0`; the last
pressure node fixes the gauge and its continuity row is dropped. `R` is a
local radius, so the operator is a frozen-radius approximation rather than a
global finite-radius stability problem.
"""
struct LSTCoefficients
    L0::Matrix{ComplexF64}       # independent of α, β, ω
    Lalpha::Matrix{ComplexF64}
    Lalpha2::Matrix{ComplexF64}
    Lbeta::Matrix{ComplexF64}
    Lbeta2::Matrix{ComplexF64}
    mass::Matrix{ComplexF64}     # coefficient of iω on the EVP right side
    keep::Vector{Int}
    n::Int
end

block(k, n) = ((k-1)*n+1):(k*n)

function _put!(matrix, row, col, values)
    matrix[row, col] .= values
end

function _coefficients(fields::AbstractMatrix, D::AbstractMatrix,
                       re_h::Real, pr::Real, R::Real)
    n = size(D, 1)
    size(D) == (n, n) || throw(DimensionMismatch("D must be square"))
    size(fields) == (4, n) || throw(DimensionMismatch("fields must be 4 × n (H,F,G,T)"))
    re_h > 0 && pr > 0 && R > 0 || throw(ArgumentError("Re_h, Pr and R must be positive"))
    H, F, G, T = (vec(fields[k, :]) for k in 1:4)
    all(isfinite, fields) && all(isfinite, D) || throw(ArgumentError("base flow must be finite"))
    sr = sqrt(re_h)
    I_n = Matrix{Float64}(I, n, n)
    D2 = D * D
    hD = Diagonal(H) * D / (R * sr)
    visc = -D2 / (R * re_h)

    L0 = zeros(ComplexF64, 5n, 5n)
    La = similar(L0); fill!(La, 0)
    La2 = similar(L0); fill!(La2, 0)
    Lb = similar(L0); fill!(Lb, 0)
    Lb2 = similar(L0); fill!(Lb2, 0)
    M = similar(L0); fill!(M, 0)
    u, v, w, theta, p = (block(k, n) for k in 1:5)

    # Momentum: advection, shear/curvature, pressure, and viscosity.
    for row in (u, v, w)
        _put!(L0, row, row, hD + visc + Diagonal(F ./ R))
        _put!(La, row, row, im * Diagonal(F))
        # The legacy local operator uses the frame rotating with the rotor.
        # G is stored as absolute swirl, so only advection uses G-1.
        _put!(Lb, row, row, im * Diagonal(G .- 1))
        _put!(La2, row, row, I_n / R)
        _put!(Lb2, row, row, I_n / R)
        _put!(M, row, row, im * I_n)
    end
    # The axial base shear is H'/(R√Re_h), whereas radial and azimuthal
    # curvature each contribute F/R. Keep the old four-field convention.
    _put!(L0, w, w, hD + visc + Diagonal((D * H) ./ (R * sr)))
    _put!(L0, u, v, -2 .* Diagonal(G) / R)
    _put!(L0, v, u,  2 .* Diagonal(G) / R)
    _put!(L0, u, w, Diagonal(D * F) / sr)
    _put!(L0, v, w, Diagonal(D * G) / sr)
    _put!(La, u, p, im * I_n)
    _put!(Lb, v, p, im * I_n)
    _put!(L0, w, p, D / sr)

    # Traditional fixed-reference centrifugal closure.  The prescribed
    # centrifugal source is linearized as -theta/R in the momentum residual.
    # There are no azimuthal or axial temperature columns in this closure.
    _put!(L0, u, theta, -I_n / R)

    # Temperature is passive when T'=0, but couples back to radial momentum.
    _put!(L0, theta, theta, hD - D2 / (R * re_h * pr))
    _put!(La, theta, theta, im * Diagonal(F))
    _put!(Lb, theta, theta, im * Diagonal(G .- 1))
    _put!(La2, theta, theta, I_n / (R * pr))
    _put!(Lb2, theta, theta, I_n / (R * pr))
    _put!(M, theta, theta, im * I_n)
    _put!(L0, theta, w, Diagonal(D * T) / sr)

    # Frozen-radius incompressibility: iαu+iβv+w_z/√Re_h+u/R=0.
    _put!(L0, p, u, I_n / R)
    _put!(L0, p, w, D / sr)
    _put!(La, p, u, im * I_n)
    _put!(Lb, p, v, im * I_n)

    removed = [first(u), last(u), first(v), last(v),
               first(w), last(w), first(theta), last(theta), last(p)]
    keep = setdiff(collect(1:5n), removed)
    return LSTCoefficients(L0, La, La2, Lb, Lb2, M, keep, n)
end

"""Assemble coefficients from a `ThermalProfile` without resolving its base flow."""
lst_coefficients(profile; radius::Real) =
    _coefficients(profile.fields, profile.D, profile.config.re_h,
                  profile.config.pr, radius)

"""Accept an existing isothermal flow for the independent legacy regression."""
lst_coefficients(flow, re_h::Real, pr::Real; radius::Real) =
    _coefficients(vcat(reshape(vec(flow.H), 1, :),
                       reshape(vec(flow.F), 1, :),
                       reshape(vec(flow.G), 1, :),
                       ones(1, length(flow.F))),
                  flow.D, re_h, pr, radius)

"""Return `(L0,L1,L2)` for `(L0+αL1+α²L2)q=0` at fixed real/complex ω."""
function spatial_matrices(cof::LSTCoefficients; beta::Real, omega::Number)
    k = cof.keep
    L0 = cof.L0 + beta * cof.Lbeta + beta^2 * cof.Lbeta2 - omega * cof.mass
    return L0[k,k], cof.Lalpha[k,k], cof.Lalpha2[k,k]
end

"""Return `(A,B)` for `Aq=ωBq`; `Im(ω)>0` means temporal growth."""
function temporal_matrices(cof::LSTCoefficients; alpha::Number, beta::Real)
    k = cof.keep
    A = cof.L0 + alpha * cof.Lalpha + alpha^2 * cof.Lalpha2 +
        beta * cof.Lbeta + beta^2 * cof.Lbeta2
    return A[k,k], cof.mass[k,k]
end

"""Restore zero wall disturbances and the fixed pressure gauge."""
function expand_mode(cof::LSTCoefficients, reduced::AbstractVector)
    length(reduced) == length(cof.keep) || throw(DimensionMismatch("reduced mode length"))
    full = zeros(ComplexF64, 5cof.n)
    full[cof.keep] .= reduced
    return (; u=full[block(1,cof.n)], v=full[block(2,cof.n)],
            w=full[block(3,cof.n)], theta=full[block(4,cof.n)],
            p=full[block(5,cof.n)])
end

end

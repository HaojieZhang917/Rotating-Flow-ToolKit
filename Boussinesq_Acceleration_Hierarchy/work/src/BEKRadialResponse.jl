module BEKRadialResponse

using LinearAlgebra
using ..BEKConsistent

export radial_response_pencil, radial_diffusion_operators

"""
Return the collocation pencil (J, Mr) for steady radial Mellin modes,
so that (J + lambda*Mr)*mode = 0.

The boundary rows of Mr vanish. J is the fixed-(Ro,Tw) similarity
Jacobian of the fully self-consistent four-field residual.
"""
function radial_response_pencil(solution::BEKConsistent.ConsistentSolution)
    op = solution.operators
    state = BEKConsistent.pack(solution.fields)
    _, J = BEKConsistent._residual_jacobian(
        state, solution.Ro, solution.Tw, solution.gamma, solution.Pr, op)

    n = length(op.x)
    F = view(state, n+1:2n)
    T = view(state, 3n+1:4n)
    Ro = solution.Ro
    chi = 1 .- solution.gamma .* (T .- 1)

    Mr = zeros(4n, 4n)
    rH, rF = 1:n, n+1:2n
    rG, rT = 2n+1:3n, 3n+1:4n
    In = Matrix{Float64}(I, n, n)
    Mr[rH, rF] .= In
    Mr[rF, rF] .= Diagonal(chi .* Ro .* F)
    Mr[rG, rG] .= Diagonal(chi .* Ro .* F)
    Mr[rT, rT] .= Diagonal(solution.Pr .* Ro .* F)

    for row in (first(rH), last(rH), first(rF), last(rF),
                first(rG), last(rG), first(rT), last(rT))
        Mr[row, :] .= 0
    end
    J, Mr
end

"""
Return the dimensionless frozen-radius radial-diffusion operators `(D2,D1)`
defined by

    J + lambda*Mr + kappa*(lambda^2*D2 + lambda*D1),

where `kappa=(ell/R0)^2`.  The cylindrical vector Laplacian gives
`F_xixi+2F_xi` and `G_xixi+2G_xi`; the scalar thermal Laplacian gives
`T_xixi`.  There is no H block because the H row in the reduced four-field
BEK system is continuity, not axial momentum.
"""
function radial_diffusion_operators(
        solution::BEKConsistent.ConsistentSolution)
    n = length(solution.operators.x)
    D2 = zeros(4n, 4n)
    D1 = zeros(4n, 4n)
    rH, rF = 1:n, n+1:2n
    rG, rT = 2n+1:3n, 3n+1:4n
    In = Matrix{Float64}(I, n, n)
    D2[rF,rF] .= In
    D2[rG,rG] .= In
    D2[rT,rT] .= In
    D1[rF,rF] .= 2In
    D1[rG,rG] .= 2In

    for row in (first(rH), last(rH), first(rF), last(rF),
                first(rG), last(rG), first(rT), last(rT))
        D2[row,:] .= 0
        D1[row,:] .= 0
    end
    D2, D1
end

end

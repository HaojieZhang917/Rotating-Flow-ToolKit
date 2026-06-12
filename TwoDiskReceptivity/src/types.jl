# Types and data structures for TwoDiskReceptivity

"""
    FlowParameters

Holds all physical and numerical parameters for the receptivity analysis.
"""
struct FlowParameters
    Res :: Int64       # Reynolds number based on disk spacing
    Ro  :: Float64     # Rotation ratio (disk angular velocity ratio)
    Ts  :: Float64     # Suction/injection parameter
    R   :: Float64     # Reynolds number
    n   :: Float64     # Azimuthal wavenumber (R * β)
    OMEGA :: Float64   # Frequency (R * ω)
    N_cheb :: Int64    # Number of Chebyshev collocation points
    mode :: Int64      # Base flow mode (1, 2, or 3)
end

"""
    ChebyshevGrid

Chebyshev collocation grid and differentiation matrices.
"""
struct ChebyshevGrid
    z  :: Vector{Float64}      # Collocation points
    D  :: Matrix{Float64}      # First derivative matrix
    D2 :: Matrix{Float64}      # Second derivative matrix
end

"""
    BaseFlowProfile

Base flow velocity profiles (F, G, H) on the Chebyshev grid.
F: scaled radial velocity
G: scaled azimuthal velocity (-1 for Bödewadt)
H: scaled wall-normal velocity
"""
struct BaseFlowProfile
    F :: Matrix{Float64}
    G :: Matrix{Float64}
    H :: Matrix{Float64}
end

"""
    CoeffMatrix

Coefficient matrices for the linearised stability equations.
These 4×4 block matrices (each block is N×N) form the building blocks
of the quadratic eigenvalue problem.
"""
struct CoeffMatrix
    Ta   :: Matrix{ComplexF64}
    A    :: Matrix{ComplexF64}
    B    :: Matrix{ComplexF64}
    C    :: Matrix{ComplexF64}
    dC   :: Matrix{ComplexF64}
    D1   :: Matrix{ComplexF64}
    Vxx  :: Matrix{ComplexF64}
    Vyy  :: Matrix{ComplexF64}
    Vzz  :: Matrix{ComplexF64}
    dVzz :: Matrix{ComplexF64}
    d2Vzz :: Matrix{ComplexF64}
    Vxy  :: Matrix{ComplexF64}
    Vxz  :: Matrix{ComplexF64}
    dVxz :: Matrix{ComplexF64}
    Vyz  :: Matrix{ComplexF64}
    dVyz :: Matrix{ComplexF64}
end

"""
    PEPMatrices

Coefficient matrices of the quadratic eigenvalue problem:
    L2 * α^2 + L1 * α + L0 = 0
where α is the complex streamwise wavenumber.
"""
struct PEPMatrices
    L0 :: Matrix{ComplexF64}
    L1 :: Matrix{ComplexF64}
    L2 :: Matrix{ComplexF64}
end

"""
    ReceptivityResult

Results of a receptivity computation.

# Fields
- `eigval_direct`: Direct eigenvalue (complex α)
- `eigval_adjoint`: Adjoint eigenvalue (complex α)
- `eigvec_direct`: Direct eigenvector
- `eigvec_adjoint`: Adjoint eigenvector
- `Cr`: Receptivity coefficient
- `BC`: Boundary condition contribution
- `Q`: Coupling coefficient (inner product)
"""
struct ReceptivityResult
    eigval_direct  :: ComplexF64
    eigval_adjoint :: ComplexF64
    eigvec_direct  :: Vector{ComplexF64}
    eigvec_adjoint :: Vector{ComplexF64}
    Cr :: Float64
    BC :: ComplexF64
    Q  :: ComplexF64
    velocity_direct  :: NTuple{4, Vector{ComplexF64}}
    velocity_adjoint :: NTuple{4, Vector{ComplexF64}}
end

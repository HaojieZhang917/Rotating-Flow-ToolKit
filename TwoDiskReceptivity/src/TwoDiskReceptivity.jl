module TwoDiskReceptivity

# Re-export key types and functions
export FlowParameters, ChebyshevGrid, BaseFlowProfile, CoeffMatrix, PEPMatrices, ReceptivityResult
export compute_baseflow, chebyshev_grid, interpolate_to_cheb, solve_baseflow
export assemble_coeff_matrix_BEK1, assemble_coeff_matrix_BEK2
export assemble_direct_matrices, assemble_adjoint_matrices
export apply_boundary_conditions!, reconstruct_eigenvector, solve_eigenvalue_problem, solve_stability_eigenproblem
export normalise_eigenvectors!, chebyshev_quadrature_weights
export compute_coupling_coefficient, compute_receptivity_coefficient, solve_receptivity

using LinearAlgebra

include("types.jl")
include("baseflow.jl")
include("stability.jl")
include("eigenvalue.jl")
include("receptivity.jl")

end # module TwoDiskReceptivity

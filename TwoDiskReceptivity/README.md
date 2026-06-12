# TwoDiskReceptivity.jl

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**TwoDiskReceptivity** is a Julia package for computing the **receptivity of compressible two-disk (rotor-stator) cavity flow**. It solves the Bödewadt–von Kármán base flow, performs direct and adjoint linear stability analysis, and computes the receptivity coefficient quantifying how wall roughness excites instability modes.

---

## Physical Background

Receptivity analysis describes how external disturbances (e.g., wall roughness, free-stream turbulence) are converted into instability waves in a boundary layer. In a two-disk configuration, the rotating disks generate a self-similar Bödewadt flow. The receptivity coefficient $C_r$ connects the wall roughness amplitude to the resulting instability mode amplitude via an adjoint-based approach.

The analysis involves three steps:

1. **Base flow**: Solve the Bödewadt–von Kármán similarity ODEs using SciPy's BVP solver.
2. **Stability**: Assemble the linearised Navier–Stokes equations into a quadratic eigenvalue problem (PEP) for the complex streamwise wavenumber $\alpha$, and solve both direct and adjoint problems.
3. **Receptivity**: Compute the coupling coefficient $Q$ (inner product of direct and adjoint modes) and the receptivity coefficient $C_r$ using the Choudhari–Street formalism.

---

## Directory Structure

```
TwoDiskReceptivity/
├── src/
│   ├── TwoDiskReceptivity.jl    # Main module
│   ├── types.jl                 # Data structures
│   ├── baseflow.jl              # Base flow solver (SciPy BVP + B-spline interpolation)
│   ├── stability.jl             # Coefficient matrix assembly (BEK1 & BEK2)
│   ├── eigenvalue.jl            # Eigenvalue problem & IAR solver
│   └── receptivity.jl           # Receptivity analysis (Q, Cr)
├── test/
│   └── runtests.jl              # Unit tests
├── examples/
│   └── example.jl               # Example usage
├── Project.toml
├── LICENSE
└── README.md
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| [BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl) | B-spline interpolation (SciPy → Chebyshev grid) |
| [PyCall.jl](https://github.com/JuliaPy/PyCall.jl) | Calling Python's SciPy `solve_bvp` |
| [NonlinearEigenproblems.jl](https://github.com/nep-pack/NonlinearEigenproblems.jl) | IAR solver for the PEP |
| [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) | Linear algebra |
| [Plots.jl](https://github.com/JuliaPlots/Plots.jl) | Plotting (optional) |
| [DelimitedFiles](https://docs.julialang.org/en/v1/stdlib/DelimitedFiles/) | I/O utilities |

**Python dependencies** (installed via PyCall / Conda): `numpy`, `scipy`

---

## Installation

```bash
git clone https://github.com/<HaojieZhang917>/TwoDiskReceptivity.jl.git
cd TwoDiskReceptivity.jl
julia --project=.
```

Then in the Julia REPL:
```julia
] instantiate
```

---

## Quick Start

```julia
using TwoDiskReceptivity

# Define parameters
params = FlowParameters(
    1000,     # Res: disk-spacing Reynolds number
    -1.0,     # Ro: rotation ratio
    0.0,      # Ts: suction/injection
    500.0,    # R: Reynolds number
    30.0,     # n: azimuthal wavenumber (R * β)
    8.0,      # OMEGA: frequency (R * ω)
    129,      # N_cheb: Chebyshev collocation points
    1         # mode: base flow type
)

# Solve base flow
baseflow, grid = solve_baseflow(params)

# Run receptivity analysis
hr = 1.0 / 1000   # roughness height
ls = 0.5           # roughness length scale
sigma = 0.4 + 0.0im

result = solve_receptivity(params, baseflow, grid, sigma, hr, ls)

println("Receptivity coefficient Cr = ", result.Cr)
println("Direct eigenvalue α = ", result.eigval_direct)
```

---

## API Reference

### Core Functions

| Function | Description |
|----------|-------------|
| `solve_baseflow(params)` | Compute base flow + Chebyshev grid |
| `solve_receptivity(params, baseflow, grid, sigma, hr, ls)` | Full receptivity pipeline → `ReceptivityResult` |
| `assemble_coeff_matrix_BEK1(...)` | Assemble coefficient matrices (mode 1 scaling) |
| `assemble_coeff_matrix_BEK2(...)` | Assemble coefficient matrices (mode 2/3 scaling) |
| `assemble_direct_matrices(cof, ...)` | Build PEP matrices L0, L1, L2 |
| `assemble_adjoint_matrices(cof, ...)` | Build adjoint PEP matrices A0, A1, A2 |
| `solve_eigenvalue_problem(L0, L1, L2, sigma, nev)` | Solve PEP via IAR |
| `normalise_eigenvectors!(...)` | Normalise direct & adjoint eigenvectors |
| `compute_coupling_coefficient(...)` | Compute coupling coefficient Q |
| `compute_receptivity_coefficient(...)` | Compute receptivity coefficient Cr |
| `chebyshev_quadrature_weights(N)` | Clenshaw-Curtis quadrature weights |

### Key Types

| Type | Fields |
|------|--------|
| `FlowParameters` | `Res, Ro, Ts, R, n, OMEGA, N_cheb, mode` |
| `ChebyshevGrid` | `z, D, D2` |
| `BaseFlowProfile` | `F, G, H` |
| `CoeffMatrix` | `Ta, A, B, C, dC, D1, Vxx, Vyy, Vzz, dVzz, d2Vzz, ...` |
| `ReceptivityResult` | `eigval_direct, eigval_adjoint, Cr, BC, Q, velocity_direct, velocity_adjoint` |

---

## Numerical Methods

1. **Base flow**: The ODE system is solved via Python SciPy's `solve_bvp`, then interpolated onto a Chebyshev collocation grid using B-splines (BSplineKit.jl).

2. **Stability analysis**: The linearised perturbation equations are discretised on the Chebyshev grid. Homogeneous Dirichlet boundary conditions are enforced by removing boundary rows/columns, yielding a **quadratic eigenvalue problem**:
   $$L_2 \alpha^2 + L_1 \alpha + L_0 = 0$$

3. **Eigenvalue solver**: The PEP is solved using the **Infinite Arnoldi Method (IAR)** with shift-and-invert (NonlinearEigenproblems.jl).

4. **Receptivity coefficient**: Using the adjoint formalism, the coupling coefficient and receptivity coefficient are computed via inner products weighted by Clenshaw-Curtis quadrature.

---

## Base Flow Modes

| mode | Description | Scaling |
|------|-------------|---------|
| 1 | Single-disk Bödewadt | Wall-normal scaled by $1/\sqrt{Re_s}$ |
| 2 | Two-disk cavity (simplified) | Fixed grid [0, 30] |
| 3 | Two-disk cavity (with κ) | Fixed grid with rotation ratio κ |

---

## Running Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## Author

- **Haojie Zhang** — [hj_zhang@tju.edu.cn](mailto:hj_zhang@tju.edu.cn)

## Citation

If you use this code in your research, please cite or acknowledge the author.

## License

This project is open-source under the MIT License. See the [LICENSE](LICENSE) file for details.

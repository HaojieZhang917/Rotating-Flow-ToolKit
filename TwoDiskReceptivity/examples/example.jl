# Example: Receptivity analysis of two-disk Bödewadt flow

using TwoDiskReceptivity
using DelimitedFiles

# ===== Physical parameters =====
Res   = 1000        # Disk-spacing Reynolds number
Ro    = -1.0        # Rotation ratio (counter-rotating)
Ts    = 0.0         # No suction/injection
R     = 500         # Reynolds number
n     = 30          # Azimuthal wavenumber
OMEGA = 8.0         # Frequency
N_cheb = 129        # Chebyshev collocation points
mode  = 1           # Base flow mode

# ===== Receptivity parameters =====
hr = 1.0 / Res      # Roughness height
ls = 0.5            # Roughness length scale
sigma = 0.4 + 0.0im # Shift-invert target

# ===== Create parameter struct =====
params = FlowParameters(Res, Ro, Ts, R, n, OMEGA, N_cheb, mode)

# ===== Solve base flow =====
println("Computing base flow...")
baseflow, grid = solve_baseflow(params)
println("  Done. Grid points: $(length(grid.z))")

# ===== Solve receptivity problem =====
println("Computing receptivity coefficient...")
result = solve_receptivity(params, baseflow, grid, sigma, hr, ls)

# ===== Output results =====
println("\n========== Results ==========")
println("Direct eigenvalue α  = $(result.eigval_direct)")
println("Adjoint eigenvalue α = $(result.eigval_adjoint)")
println("Coupling coefficient |Q| = $(abs(result.Q))")
println("Boundary term |BC|     = $(abs(result.BC))")
println("Receptivity coefficient Cr = $(result.Cr)")
println("==============================")

# ===== Save eigenfunctions =====
vel_dir = result.velocity_direct
vel_adj = result.velocity_adjoint
z = grid.z

writedlm("eigenfunctions.dat",
    [z abs.(vel_dir[1]) abs.(vel_dir[2]) abs.(vel_dir[3]) abs.(vel_dir[4])
        abs.(vel_adj[1]) abs.(vel_adj[2]) abs.(vel_adj[3]) abs.(vel_adj[4])])
println("\nEigenfunctions saved to eigenfunctions.dat")

using LinearAlgebra
using Printf
using DelimitedFiles

# ---------------- User parameters ----------------
const N = 60
const Reh = 2000.0
const CQ = 0.0
const r0 = 0.05
const Ts = 0.0
const output_file = joinpath(@__DIR__, "throughflow_initial_profile.csv")

const newton_tol = 1e-10
const newton_maxit = 50

"""
    cheb_lobatto()

Chebyshev-Gauss-Lobatto points on z in [0, 1], first derivative matrix D,
second derivative matrix D2, and Clenshaw-Curtis quadrature weights w.

The points are returned in increasing order: z[1] = 0, z[end] = 1.
"""
function cheb_lobatto()
    N < 2 && error("N must be at least 2")

    x = [cos(pi * j / N) for j in 0:N]
    c = [j == 0 || j == N ? 2.0 : 1.0 for j in 0:N]
    c = c .* [(-1.0)^j for j in 0:N]

    D = zeros(N + 1, N + 1)
    for i in 1:N + 1
        for j in 1:N + 1
            if i != j
                D[i, j] = (c[i] / c[j]) / (x[i] - x[j])
            end
        end
    end
    for i in 1:N + 1
        D[i, i] = -sum(D[i, :])
    end

    # Reverse to increasing z and map x in [-1, 1] to z in [0, 1].
    perm = collect(N + 1:-1:1)
    x = x[perm]
    D = D[perm, perm]
    z = (x .+ 1.0) ./ 2.0

    # dx/dz = 2, so d/dz = 2 d/dx.
    Dz = 2.0 .* D
    D2z = Dz * Dz

    deriv_error = norm(Dz * z .- 1.0, Inf)
    deriv_error > 1e-8 && error("Chebyshev D matrix failed Dz*z=1 check")

    w = clenshaw_curtis_weights() ./ 2.0
    return z, Dz, D2z, w[perm]
end

"""
    clenshaw_curtis_weights()

Full-interval Clenshaw-Curtis quadrature weights on x in [-1, 1].
"""
function clenshaw_curtis_weights()
    w = zeros(N + 1)
    for j in 0:N
        s = 0.0
        for k in 1:floor(Int, N / 2)
            term = 2.0 / (1.0 - (2k)^2) * cos(2k * j * pi / N)
            s += (2k == N) ? 0.5 * term : term
        end
        cj = (j == 0 || j == N) ? 1.0 : 2.0
        w[j + 1] = cj / N * (1.0 + s)
    end
    return w
end

"""
    residual_base(x, D)

Zero-throughflow rotor-stator base flow in the same first-order form as the
user's PyCall `solve_bvp` code.

Unknown vector:

    H, F, Fz, Fzz, G, Gz

ODE:

    H_z    = -2 sqrt(Reh) F
    F_z    = Fz
    Fz_z   = Fzz
    Fzz_z  = Reh * (H Fzz/sqrt(Reh) + H_z Fz/sqrt(Reh) - 2G Gz + 2F Fz)
    G_z    = Gz
    Gz_z   = Reh * (H Gz/sqrt(Reh) + 2F G)

Boundary conditions:

    z = 0 rotor: H=-Ts, F=0, G=1
    z = 1 stator: H=0,   F=0, G=0
"""
function residual_base(x::Vector{Float64}, D)
    n = length(D[:, 1])
    H = x[1:n]
    F = x[n + 1:2n]
    Fz_var = x[2n + 1:3n]
    Fzz_var = x[3n + 1:4n]
    G = x[4n + 1:5n]
    Gz_var = x[5n + 1:6n]

    sqr = sqrt(Reh)
    Hz = -2.0 .* sqr .* F

    R = zeros(6n)
    R[1:n] = D * H .- Hz
    R[n + 1:2n] = D * F .- Fz_var
    R[2n + 1:3n] = D * Fz_var .- Fzz_var
    R[3n + 1:4n] = D * Fzz_var .-
        Reh .* ((H .* Fzz_var .+ Fz_var .* Hz) ./ sqr .- 2.0 .* G .* Gz_var .+ 2.0 .* F .* Fz_var)
    R[4n + 1:5n] = D * G .- Gz_var
    R[5n + 1:6n] = D * Gz_var .-
        Reh .* ((H .* Gz_var) ./ sqr .+ 2.0 .* F .* G)

    # Boundary conditions replace selected tau rows.
    R[1] = H[1] + Ts
    R[n] = H[end]
    R[n + 1] = F[1]
    R[2n] = F[end]
    R[4n + 1] = G[1] - 1.0
    R[5n] = G[end]
    return R
end

function finite_difference_jacobian(f, x, eps_scale)
    fx = f(x)
    J = zeros(length(fx), length(x))
    for j in eachindex(x)
        h = eps_scale * max(1.0, abs(x[j]))
        xp = copy(x)
        xp[j] += h
        J[:, j] = (f(xp) - fx) ./ h
    end
    return J, fx
end

function solve_base_flow()
    z, D, D2, w = cheb_lobatto()
    n = N + 1

    # Smooth initial guess satisfying wall values.
    H = -Ts .* (1.0 .- z)
    F = zeros(n)
    Fz_var = zeros(n)
    Fzz_var = zeros(n)
    G = 1.0 .- z
    Gz_var = -ones(n)
    x = vcat(H, F, Fz_var, Fzz_var, G, Gz_var)

    f = y -> residual_base(y, D)
    final_res = Inf

    for it in 1:newton_maxit
        J, R = finite_difference_jacobian(f, x, sqrt(eps(Float64)))
        res = norm(R, Inf)
        final_res = res
        @printf("base Newton %2d: ||R||_inf = %.4e\n", it, res)
        res < newton_tol && break

        dx = -(J \ R)
        alpha = 1.0
        old = res
        while alpha > 1e-4
            trial = x .+ alpha .* dx
            new = norm(f(trial), Inf)
            if new < old
                x = trial
                break
            end
            alpha *= 0.5
        end
        alpha <= 1e-4 && error("Base-flow Newton line search failed")
    end

    final_res = norm(f(x), Inf)
    final_res > newton_tol && error("Base-flow Newton did not converge. Final ||R||_inf = $(final_res)")

    H = x[1:n]
    F = x[n + 1:2n]
    Fz_var = x[2n + 1:3n]
    Fzz_var = x[3n + 1:4n]
    G = x[4n + 1:5n]

    K0_values = Fzz_var ./ Reh .- Fz_var .* H ./ sqrt(Reh) .+ (G .^ 2 .- F .^ 2)
    K0 = -sum(w .* K0_values)
    return z, D, D2, w, F, G, H, K0, final_res
end

"""
    solve_throughflow_correction(F, G, H, D, D2, w)

Solve the first center-throughflow correction:

    F1_zz/Reh - H F1_z/sqrt(Reh) + 2G G1 + K = 0
    G1_zz/Reh - H G1_z/sqrt(Reh) - 2G F1 = 0

with no-slip correction at both walls:

    F1(0)=F1(1)=G1(0)=G1(1)=0

and radial flux constraint:

    integral_0^1 F1 dz = CQ

Here K is the constant first-order radial pressure-gradient correction.
For this expansion H1(z) = 0.
"""
function solve_throughflow_correction(F, G, H, D, D2, w)
    n = length(F)
    sqr = sqrt(Reh)

    A = zeros(2n + 1, 2n + 1)
    b = zeros(2n + 1)

    # Unknown vector y = [F1; G1; K].
    for i in 1:n
        # F1 equation
        A[i, 1:n] .= D2[i, :] ./ Reh .- H[i] .* D[i, :] ./ sqr
        A[i, n + i] += 2.0 * G[i]
        A[i, 2n + 1] = 1.0

        # G1 equation
        row = n + i
        A[row, n + 1:2n] .= D2[i, :] ./ Reh .- H[i] .* D[i, :] ./ sqr
        A[row, i] += -2.0 * G[i]
    end

    # Boundary rows.
    A[1, :] .= 0.0
    A[1, 1] = 1.0
    b[1] = 0.0

    A[n, :] .= 0.0
    A[n, n] = 1.0
    b[n] = 0.0

    A[n + 1, :] .= 0.0
    A[n + 1, n + 1] = 1.0
    b[n + 1] = 0.0

    A[2n, :] .= 0.0
    A[2n, 2n] = 1.0
    b[2n] = 0.0

    # Integral flux constraint.
    A[2n + 1, :] .= 0.0
    A[2n + 1, 1:n] .= w
    b[2n + 1] = CQ

    y = A \ b
    linear_residual = norm(A * y - b, Inf)
    matrix_condition = cond(A)
    F1 = y[1:n]
    G1 = y[n + 1:2n]
    K = y[2n + 1]
    H1 = zeros(n)
    return F1, G1, H1, K, linear_residual, matrix_condition
end

function reconstruct_pressure(z, D, D2, H)
    sqr = sqrt(Reh)
    Hz = D * H
    Hzz = D2 * H
    Pz = H .* Hz .- Hzz ./ sqr

    # Simple cumulative trapezoid for a plotted pressure profile, P(0)=0.
    P = zeros(length(z))
    for i in 2:length(z)
        P[i] = P[i - 1] + 0.5 * (Pz[i] + Pz[i - 1]) * (z[i] - z[i - 1])
    end
    return P, Pz
end

function compute_profile()
    r0 <= 0 && error("r0 must be positive. A finite center flux is singular at exactly r = 0.")

    z, D, D2, w, F, G, H, K0, base_residual = solve_base_flow()
    F1, G1, H1, K, linear_residual, matrix_condition = solve_throughflow_correction(F, G, H, D, D2, w)
    P, Pz = reconstruct_pressure(z, D, D2, H)

    # Physical/asymptotic profile near the center at r = r0.
    U = r0 .* F .+ F1 ./ r0
    V = r0 .* G .+ G1 ./ r0
    W = H .+ H1 ./ r0^2

    data = hcat(z, F, G, H, F1, G1, H1, U, V, W, P, Pz)
    header = "z,F0,G0,H0,F1,G1,H1,U_at_r0,V_at_r0,W_at_r0,P0,Pz0"
    open(output_file, "w") do io
        println(io, header)
        writedlm(io, data, ',')
    end

    @printf("\nDone.\n")
    @printf("Reh = %.6g, CQ = %.6g, r0 = %.6g\n", Reh, CQ, r0)
    @printf("Ts = %.6e\n", Ts)
    @printf("Solved zero-order radial pressure constant K0 = %.6e\n", K0)
    @printf("Zero-order nonlinear residual = %.6e\n", base_residual)
    @printf("Integral(F0) = %.6e\n", dot(w, F))
    @printf("Integral(F1) = %.6e  target CQ = %.6e\n", dot(w, F1), CQ)
    @printf("First-order linear residual = %.6e\n", linear_residual)
    @printf("First-order matrix condition = %.6e\n", matrix_condition)
    @printf("First-order pressure-gradient correction K = %.6e\n", K)
    @printf("Output written to %s\n", output_file)

    return (; z, F, G, H, F1, G1, H1, U, V, W, P, Pz, K0, K, w, base_residual,
            linear_residual, matrix_condition)
end

if abspath(PROGRAM_FILE) == @__FILE__
    compute_profile()
end

"""
    ChebNewton.jl — Chebyshev 配置法 + 解析 Jacobian + Newton 迭代
    求解动坐标系 von Kármán 旋转圆盘不可压缩 Boussinesq 基本流

    用法:
        include("ChebNewton.jl")
        Y, x = solve_baseflow(N=40, Tw=1.0)
"""
module ChebNewton
using LinearAlgebra, Printf

export solve_baseflow, cheb, map_interval

# ═══════════════════════════════════════════════════════
#  Chebyshev 微分矩阵 ([-1,1])
# ═══════════════════════════════════════════════════════
function cheb(N::Int)
    N == 0 && return [0.0], [0.0;;]
    xi = cos.(π * (0:N) / N)
    c = ones(N+1); c[1] = 2; c[end] = 2
    c .*= (-1) .^ (0:N)
    X = repeat(xi, 1, N+1)
    dX = X .- X' + I(N+1)
    D = (c ./ c') ./ dX
    D .-= diagm(vec(sum(D, dims=2)))
    return xi, D
end

# ═══════════════════════════════════════════════════════
#  线性映射 [-1,1] → [a,b]
# ═══════════════════════════════════════════════════════
function map_interval(N::Int, a::Real, b::Real)
    xi, Dc = cheb(N)
    x = @. a + (b - a) * (1 - xi) / 2
    D = -(2.0 / (b - a)) * Dc
    return x, D
end

# ═══════════════════════════════════════════════════════
#  残差 R(Y) = D⊗I·Y - f(Y) + BC
#  Y: (N+1,7) → R: (N+1,7)
# ═══════════════════════════════════════════════════════
function residual!(R, Y, D, Pr, Tw)
    M = size(Y, 1)
    DY = D * Y                     # 导数矩阵

    for i in 1:M
        y = view(Y, i, :)
        # ODE 残差
        R[i,1] = DY[i,1] + 2*y[3]                                     # H' + 2F
        R[i,2] = DY[i,2] - y[3]^2 - y[1]*y[2] + (y[5]-1)^2 + y[7] - 1 # F''
        R[i,3] = DY[i,3] - y[2]                                        # F'
        R[i,4] = DY[i,4] - 2*y[3]*y[5] - y[1]*y[4] + 2*y[3]           # G''
        R[i,5] = DY[i,5] - y[4]                                        # G'
        R[i,6] = DY[i,6] - Pr * y[1] * y[6]                            # T''
        R[i,7] = DY[i,7] - y[6]                                        # T'
    end

    # BCs: 壁面 z=0
    R[1,1] = Y[1,1]          # H(0) = 0
    R[1,3] = Y[1,3]          # F(0) = 0
    R[1,5] = Y[1,5]          # G(0) = 0
    R[1,7] = Y[1,7] - Tw     # T(0) = Tw

    # BCs: 远场 z=zmax
    R[M,3] = Y[M,3]          # F(∞) = 0
    R[M,5] = Y[M,5] - 1.0    # G(∞) = 1
    R[M,7] = Y[M,7] - 1.0    # T(∞) = 1
    return R
end

# ═══════════════════════════════════════════════════════
#  解析 Jacobian: J = D⊗I₇ - diag(∂f/∂Y)
#  BC 行替换为单位行
# ═══════════════════════════════════════════════════════
function build_jacobian!(J, Y, D, Pr, Tw)
    M, nv = size(Y)            # M = N+1, nv = 7
    Ntot = M * nv
    fill!(J, 0.0)

    # Part 1: D⊗I₇ (导数耦合)
    for k in 1:M, i in 1:M
        Dik = D[i,k]
        abs(Dik) < 1e-14 && continue
        for j in 1:nv
            J[(i-1)*nv+j, (k-1)*nv+j] = Dik
        end
    end

    # Part 2: -∂f/∂Y (局部非线性 Jacobian)
    for i in 1:M
        y = view(Y, i, :)
        r0 = (i-1) * nv

        # R[0] = D@H + 2F → df0/dF = -2
        J[r0+1, r0+3] -= -2.0

        # R[1] = D@F' - f1, f1 = F² + H·F' - (G-1)² - (T-1)
        J[r0+2, r0+1] -= y[2]           # -df1/dH = -F'
        J[r0+2, r0+2] -= y[1]           # -df1/dF' = -H
        J[r0+2, r0+3] -= 2*y[3]         # -df1/dF = -2F
        J[r0+2, r0+5] -= -2*(y[5]-1)    # -df1/dG = +2(G-1)
        J[r0+2, r0+7] -= -1.0           # -df1/dT = +1

        # R[2] = D@F - F' → df2/dF' = -1
        J[r0+3, r0+2] -= 1.0

        # R[3] = D@G' - f3, f3 = 2FG + H·G' - 2F
        J[r0+4, r0+1] -= y[4]           # -df3/dH = -G'
        J[r0+4, r0+3] -= 2*y[5] - 2     # -df3/dF = -(2G-2)
        J[r0+4, r0+4] -= y[1]           # -df3/dG' = -H
        J[r0+4, r0+5] -= 2*y[3]         # -df3/dG = -2F

        # R[4] = D@G - G'
        J[r0+5, r0+4] -= 1.0

        # R[5] = D@T' - Pr·H·T'
        J[r0+6, r0+1] -= Pr*y[6]        # -df5/dH = -Pr·T'
        J[r0+6, r0+6] -= Pr*y[1]        # -df5/dT' = -Pr·H

        # R[6] = D@T - T'
        J[r0+7, r0+6] -= 1.0
    end

    # Part 3: 替换 BC 行为单位矩阵
    bc_rows = [1, 3, 5, 7]           # z=0: H, F, G, T
    for j in bc_rows
        row = j
        for c in 1:Ntot
            J[row,c] = 0.0
        end
        J[row,row] = 1.0
    end
    for j in [3, 5, 7]               # z=zmax: F, G, T
        row = (M-1)*nv + j
        for c in 1:Ntot
            J[row,c] = 0.0
        end
        J[row,row] = 1.0
    end

    return J
end

# ═══════════════════════════════════════════════════════
#  Newton 迭代主函数
# ═══════════════════════════════════════════════════════
function solve_baseflow(; N::Int=40, zmax::Real=20.0, Pr::Real=0.72, Tw::Real=1.0,
                         max_iter::Int=30, rtol::Real=1e-8, verbose::Bool=true)
    # Chebyshev 网格
    x, D = map_interval(N, 0.0, zmax)
    M = N + 1
    nv = 7
    Ntot = M * nv

    # 初始猜测 (经典 von Kármán 不可压解)
    Y = zeros(M, nv)
    Y[:,1] = @. -0.8845 * (1 - exp(-0.8*x))           # H
    Y[:,2] = @.  0.5102 * (1 - 0.8*x) * exp(-0.8*x)   # F'
    Y[:,3] = @.  0.5102 * x * exp(-0.8*x)              # F
    Y[:,4] = @.  0.6159 * exp(-0.8*x)                  # G'
    Y[:,5] = @.  1.0 - exp(-0.8*x)                     # G
    Y[:,6] .= 0.0                                       # T'
    Y[:,7] .= 1.0                                       # T

    R = similar(Y)
    J = zeros(Ntot, Ntot)

    for iter in 1:max_iter
        residual!(R, Y, D, Pr, Tw)
        rnorm = norm(R)
        verbose && println("  iter $iter: ||R|| = $(@sprintf("%.2e", rnorm))")
        rnorm < rtol && break

        build_jacobian!(J, Y, D, Pr, Tw)
        dY_flat = J \ (-vec(R))
        dY = reshape(dY_flat, M, nv)

        # 阻尼线搜索
        alpha = 1.0
        Y_new = Y + alpha * dY
        R_new = similar(R)
        residual!(R_new, Y_new, D, Pr, Tw)
        r_new = norm(R_new)
        while r_new > rnorm && alpha > 1e-4
            alpha *= 0.5
            Y_new = Y + alpha * dY
            residual!(R_new, Y_new, D, Pr, Tw)
            r_new = norm(R_new)
        end
        Y = Y_new
    end

    return Y, x, D
end

end  # module

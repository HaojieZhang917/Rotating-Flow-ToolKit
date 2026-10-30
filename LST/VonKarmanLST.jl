"""
VonKarmanLST.jl
==============

自包含的 von Kármán 旋转圆盘 LST 求解工具。

本文件不 include 任何其他项目脚本。它自己完成：
1. von Kármán 基流 ODE 求解和 Hermite 插值；
2. Chebyshev 微分矩阵和有限/半无限坐标映射；
3. 时间模式矩阵算子组装；
4. 二次空间模式矩阵算子组装；
5. 边界条件处理和广义特征值求解。

时间模式：固定 alpha、beta、R，求 omega。
空间模式：固定 omega、beta、R，求 alpha。
"""
module VonKarmanLST

using LinearAlgebra

export LSTParameters, BaseFlow, VonKarmanBaseParameters,
       solve_von_karman_baseflow,
       finite_grid, affine_grid, semi_infinite_grid,
       assemble_temporal, solve_temporal,
       assemble_spatial, solve_spatial,
       SpatialOperatorBlocks, assemble_spatial_blocks,
       assemble_adjoint_spatial, solve_adjoint_spatial,
       assemble_alst, solve_alst,
       assemble_forced_spatial, solve_forced_spatial,
       wall_shear_response, solve_wall_forced_spatial,
       roughness_wall_scattering, solve_roughness_distortion,
       gaussian_roughness, roughness_fourier_spectrum,
       roughness_fourier_grid, reconstruct_roughness_field,
       reconstruct_coupled_field, nonlinear_cross_forcing,
       solve_coupled_response,
       clenshaw_curtis_weights, subinterval_cc_grid,
       vortical_polarizations, vortical_edge_states,
       response_gains, farfield_matrix, solve_farfield_modes,
       farfield_state_pencil, solve_farfield_state_modes,
       farfield_bounded_subspace, farfield_state_boundary_matrix,
       farfield_algebraic_constraints,
       farfield_mode_diagnostics,
       adjoint_overlap, compare_lsts_alst,
       solve_lsts

Base.@kwdef struct LSTParameters
    mode::Symbol = :temporal
    R::Float64 = 285.36
    alpha::Float64 = 0.38482
    beta::Float64 = 0.07759
    omega::Float64 = 0.0
    # Ro、Co 沿用 BEK 旋转项约定；当前严格数值基准为 Ro=1、Co=2
    Ro::Float64 = 1.0
    Co::Float64 = 2.0
    N::Int = 99
    temporal_zmax::Float64 = 50.0
    spatial_shift::ComplexF64 = 0.38 + 0im
    nev::Int = 5
end

struct BaseFlow
    z::Vector{Float64}
    w::Vector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
    du::Vector{Float64}
    dv::Vector{Float64}
    dw::Vector{Float64}
end

Base.@kwdef struct VonKarmanBaseParameters
    # 空间模式使用半无限映射；尾部必须足够接近远场，避免端点夹断误差
    eta_max::Float64 = 30.0
    n_points::Int = 6001
    initial_fprime::Float64 = 0.51
    initial_gprime::Float64 = -0.616
    newton_tolerance::Float64 = 1e-11
    newton_max_iterations::Int = 20
end

"""von Kármán 相似性方程的一阶形式。"""
function _von_karman_rhs(y)
    f, fp, g, gp, h = y
    return [fp, f^2 - g^2 + h*fp, gp, 2*f*g + h*gp, -2*f]
end

"""用固定步长 RK4 从 eta=0 积分到 eta_max。"""
function _integrate_von_karman(fprime0, gprime0, bp::VonKarmanBaseParameters;
                                store::Bool=false)
    n = bp.n_points
    n >= 101 || error("n_points 至少应为 101")
    eta = collect(range(0.0, bp.eta_max; length=n))
    step = eta[2] - eta[1]
    y = [0.0, fprime0, 1.0, gprime0, 0.0]
    states = store ? Matrix{Float64}(undef, 5, n) : Matrix{Float64}(undef, 0, 0)
    store && (states[:, 1] .= y)
    for j in 1:n-1
        k1 = _von_karman_rhs(y)
        k2 = _von_karman_rhs(y .+ 0.5*step.*k1)
        k3 = _von_karman_rhs(y .+ 0.5*step.*k2)
        k4 = _von_karman_rhs(y .+ step.*k3)
        y .+= (step/6.0).*(k1 .+ 2k2 .+ 2k3 .+ k4)
        all(isfinite, y) || error("von Kármán ODE 积分发散，请检查初值或 eta_max")
        store && (states[:, j+1] .= y)
    end
    return store ? (eta, states) : y
end

"""求解 von Kármán 相似性基本流，并返回随盘坐标系速度。"""
function solve_von_karman_baseflow(; params::VonKarmanBaseParameters=VonKarmanBaseParameters())
    a = params.initial_fprime
    b = params.initial_gprime
    converged = false
    for _ in 1:params.newton_max_iterations
        y = _integrate_von_karman(a, b, params)
        residual = [y[1], y[3]]
        if norm(residual, Inf) < params.newton_tolerance
            converged = true
            break
        end
        da = 1e-5 * max(1.0, abs(a))
        db = 1e-5 * max(1.0, abs(b))
        ya = _integrate_von_karman(a + da, b, params)
        yb = _integrate_von_karman(a, b + db, params)
        J = [(ya[1]-y[1])/da (yb[1]-y[1])/db;
             (ya[3]-y[3])/da (yb[3]-y[3])/db]
        delta = J \ (-residual)
        a += delta[1]
        b += delta[2]
    end
    converged || error("von Kármán 基流 Newton 射击未收敛")
    eta, states = _integrate_von_karman(a, b, params; store=true)
    f, fp, g, gp, h = states[1,:], states[2,:], states[3,:], states[4,:], states[5,:]
    return BaseFlow(eta, h, f, g .- 1.0, fp, gp, -2.0 .* f)
end

function _hermite_value(x, y, dy, q)
    q <= x[1] && return y[1]
    q >= x[end] && return y[end]
    i = min(max(searchsortedlast(x, q), 1), length(x)-1)
    h = x[i+1] - x[i]
    t = (q-x[i])/h
    t2 = t*t; t3 = t2*t
    return (2t3-3t2+1)*y[i] + (t3-2t2+t)*h*dy[i] +
           (-2t3+3t2)*y[i+1] + (t3-t2)*h*dy[i+1]
end

function _hermite_pair(x, y, dy, q)
    if q <= x[1] || q >= x[end]
        return q <= x[1] ? (y[1], 0.0) : (y[end], 0.0)
    end
    i = min(max(searchsortedlast(x, q), 1), length(x)-1)
    h = x[i+1] - x[i]
    t = (q-x[i])/h; t2=t*t; t3=t2*t
    val = (2t3-3t2+1)*y[i] + (t3-2t2+t)*h*dy[i] +
          (-2t3+3t2)*y[i+1] + (t3-t2)*h*dy[i+1]
    der = ((6t2-6t)/h)*y[i] + (3t2-4t+1)*dy[i] +
          ((-6t2+6t)/h)*y[i+1] + (3t2-2t)*dy[i+1]
    return val, der
end

function _profile_at(base::BaseFlow, zq)
    u, du = _hermite_pair(base.z, base.u, base.du, zq)
    v, dv = _hermite_pair(base.z, base.v, base.dv, zq)
    w = _hermite_value(base.z, base.w, base.dw, zq)
    return u, v, w, du, dv, -2u
end

function _chebyshev_matrix(N::Int)
    N >= 8 || error("Chebyshev 阶数 N 至少取 8")
    x = -cos.(range(0.0, pi; length=N+1))
    c = [2.0; ones(N-1); 2.0] .* (-1.0).^(0:N)
    X = repeat(reshape(x, :, 1), 1, N+1)
    dX = X - X'
    D = (c*(1.0 ./ c)') ./ (dX + Matrix{Float64}(I, N+1, N+1))
    D -= Diagonal(vec(sum(D, dims=2)))
    return x, D
end

"""有限截断映射，时间模式使用。"""
function finite_grid(N::Int, zmax::Float64)
    x, Dx = _chebyshev_matrix(N)
    z = zeros(N+1); scale = zeros(N+1)
    for j in eachindex(x)
        xi=x[j]
        den=-2xi^3+xi^2-3xi+4
        z[j]=min(max((4xi^3-2xi^2+6xi+12)/den,0.0),zmax)
        scale[j]=(2xi^3-xi^2+3xi-4)^2/(20*(6xi^2-2xi+3))
    end
    return z, Diagonal(scale)*Dx
end

"""严格映射到 [0,zmax] 的仿射 Chebyshev 网格，受迫外缘问题使用。"""
function affine_grid(N::Int,zmax::Float64)
    zmax>0 || error("zmax 必须为正")
    x,Dx=_chebyshev_matrix(N)
    z=(zmax/2.0).*(x .+ 1.0)
    D=(2.0/zmax).*Dx
    return z,D
end

"""
Clenshaw–Curtis 权重，节点顺序为 `x=-cos(j*pi/N)`。

权重对应整个 [-1,1] 区间；子区间积分必须重新生成节点和 Jacobian，
不能简单删去原网格中位于区间外的节点。
"""
function clenshaw_curtis_weights(N::Int)
    N >= 1 || error("Clenshaw–Curtis 阶数 N 必须为正")
    w=zeros(Float64,N+1)
    for j in 0:N
        s=0.0
        for k in 1:floor(Int,N/2)
            term=2.0/(1.0-(2k)^2)*cos(2k*j*pi/N)
            s += (2k==N ? 0.5*term : term)
        end
        endpoint=(j==0 || j==N) ? 1.0 : 2.0
        w[j+1]=(endpoint/N)*(1.0+s)
    end
    return w
end

"""[0,zi] 上重新生成的 Clenshaw–Curtis 节点和物理权重。"""
function subinterval_cc_grid(zi::Float64,Nq::Int=128)
    zi>0 || error("zi 必须为正")
    ξ=-cos.(range(0.0,pi;length=Nq+1))
    z=(zi/2.0).*(ξ .+ 1.0)
    w=(zi/2.0).*clenshaw_curtis_weights(Nq)
    return z,w
end

"""给定波矢的两个正交不可压涡波极化向量。"""
function vortical_polarizations(k::AbstractVector{<:Real})
    length(k)==3 || error("k 必须包含三个实波数分量")
    kk=Float64.(collect(k)); kn=norm(kk)
    kn>0 || error("波矢不能为零")
    khat=kk/kn
    ref=abs(dot(khat,[0.0,0.0,1.0]))<0.9 ? [0.0,0.0,1.0] : [1.0,0.0,0.0]
    e1=cross(khat,ref); e1/=norm(e1)
    e2=cross(khat,e1); e2/=norm(e2)
    dot(e2,ref)<0 && (e2=-e2)
    return hcat(e1,e2)
end

"""用两个极化向量构造等输入速度幅值的四变量外缘状态。"""
function vortical_edge_states(k::AbstractVector{<:Real};
                              amplitude::Real=1.0,
                              pressure::Complex=0.0+0im)
    E=vortical_polarizations(k)
    states=Vector{Vector{ComplexF64}}(undef,2)
    for j in 1:2
        states[j]=ComplexF64[amplitude*E[1,j],amplitude*E[2,j],
                             amplitude*E[3,j],pressure]
    end
    return states
end

"""
组装当前旋转坐标系 von Kármán 远场的 4×4 常系数矩阵。

令法向依赖为 `exp(im*kz*z)`，则

    M_inf(kz) = M0 + kz*M1 + kz^2*M2

满足 `M_inf(kz)*a = 0` 才是当前局部原始变量算子的一致远场模态。
"""
function farfield_matrix(base::BaseFlow;
                         p::LSTParameters=LSTParameters(mode=:spatial),
                         kz::Number=0.0+0im)
    F=ComplexF64(base.u[end]); G=ComplexF64(base.v[end]); H=ComplexF64(base.w[end])
    dF=ComplexF64(base.du[end]); dG=ComplexF64(base.dv[end]); dH=-2F
    R=p.R; Ro=p.Ro; Co=p.Co; Z4=zeros(ComplexF64,4,4)
    Ta=ComplexF64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
    A=ComplexF64[F 0 0 1; 0 F 0 0; 0 0 F 0; 1 0 0 0]
    B=ComplexF64[G/R 0 0 0; 0 G/R 0 1/R; 0 0 G/R 0; 0 1/R 0 0]
    C=ComplexF64[Ro*H/R 0 0 0; 0 Ro*H/R 0 0; 0 0 Ro*H/R 1; 0 0 1 0]
    D1=ComplexF64[Ro*F/R -(2Ro*G+Co)/R dF 0;
                  (2Ro*G+Co)/R Ro*F/R dG 0;
                  0 0 Ro*dH/R 0;
                  Ro/R 0 0 0]
    Vxx=ComplexF64[-1/R 0 0 0; 0 -1/R 0 0; 0 0 -1/R 0; 0 0 0 0]
    Vyy=ComplexF64[-1/R^3 0 0 0; 0 -1/R^3 0 0; 0 0 -1/R^3 0; 0 0 0 0]
    Vzz=ComplexF64[-1/R 0 0 0; 0 -1/R 0 0; 0 0 -1/R 0; 0 0 0 0]
    M0=D1+im*p.beta*R*B-im*p.omega*Ta-p.beta^2*R^2*Vyy+
       im*p.alpha*A-p.alpha^2*Vxx
    M1=im*C
    M2=-Vzz
    return M0+ComplexF64(kz)*M1+ComplexF64(kz)^2*M2
end

"""求解当前远场 4×4 常系数算子的所有有限法向波数模态。"""
function solve_farfield_modes(base::BaseFlow;
                              p::LSTParameters=LSTParameters(mode=:spatial),
                              max_abs_kz::Real=100.0)
    F=ComplexF64(base.u[end]); G=ComplexF64(base.v[end]); H=ComplexF64(base.w[end])
    dF=ComplexF64(base.du[end]); dG=ComplexF64(base.dv[end]); dH=-2F
    R=p.R; Ro=p.Ro; Co=p.Co; Z4=zeros(ComplexF64,4,4)
    Ta=ComplexF64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0]
    A0=ComplexF64[F 0 0 1; 0 F 0 0; 0 0 F 0; 1 0 0 0]
    B0=ComplexF64[G/R 0 0 0; 0 G/R 0 1/R; 0 0 G/R 0; 0 1/R 0 0]
    C0=ComplexF64[Ro*H/R 0 0 0; 0 Ro*H/R 0 0; 0 0 Ro*H/R 1; 0 0 1 0]
    D1=ComplexF64[Ro*F/R -(2Ro*G+Co)/R dF 0;
                  (2Ro*G+Co)/R Ro*F/R dG 0;
                  0 0 Ro*dH/R 0;
                  Ro/R 0 0 0]
    Vxx=ComplexF64[-1/R 0 0 0; 0 -1/R 0 0; 0 0 -1/R 0; 0 0 0 0]
    Vyy=ComplexF64[-1/R^3 0 0 0; 0 -1/R^3 0 0; 0 0 -1/R^3 0; 0 0 0 0]
    Vzz=ComplexF64[-1/R 0 0 0; 0 -1/R 0 0; 0 0 -1/R 0; 0 0 0 0]
    M0=D1+im*p.beta*R*B0-im*p.omega*Ta-p.beta^2*R^2*Vyy+
       im*p.alpha*A0-p.alpha^2*Vxx
    M1=im*C0; M2=-Vzz
    I4=Matrix{ComplexF64}(I,4,4)
    pencil_A=[-M1 -M0; I4 Z4]
    pencil_B=[M2 Z4; Z4 I4]
    E=eigen(pencil_A,pencil_B)
    vals=ComplexF64.(E.values); vecs=ComplexF64.(E.vectors)
    ids=findall(i->isfinite(real(vals[i])) && isfinite(imag(vals[i])) &&
                         abs(vals[i])<=max_abs_kz,eachindex(vals))
    sort!(ids,by=i->abs(vals[i]))
    return [(kz=vals[i],amplitude=vecs[5:8,i],
             residual=farfield_matrix(base;p=p,kz=vals[i])*vecs[5:8,i])
            for i in ids]
end

"""
构造远场一阶广义状态系统。

当前原始变量远场方程写成

    M0*q + M1*qz_wave + M2*qz_wave^2 = 0,

其中 `qz_wave` 表示法向波数 `kz`。令实际法向导数状态为
`r=∂q/∂z`，并采用 `q(z)=a*exp(i*kz*z)`，则 `r=i*kz*q`。完整状态
`X=[q;r]` 满足奇异广义系统

    E*X_z = A*X,

其中奇异性来自压力的代数约束。保留奇异 pencil 很重要，因为直接
消去压力会丢失两个无限特征值所代表的约束信息。
"""
function farfield_state_pencil(base::BaseFlow;
                               p::LSTParameters=LSTParameters(mode=:spatial))
    M0=farfield_matrix(base;p=p,kz=0.0+0im)
    Mp=farfield_matrix(base;p=p,kz=1.0+0im)
    Mm=farfield_matrix(base;p=p,kz=-1.0+0im)
    M1=(Mp-Mm)/2
    M2=(Mp+Mm-2M0)/2
    Z4=zeros(ComplexF64,4,4)
    I4=Matrix{ComplexF64}(I,4,4)
    E=[I4 Z4; Z4 M2]
    A=[Z4 I4; M0 -im*M1]
    return (E=E,A=A,M0=M0,M1=M1,M2=M2,params=p)
end

"""
求解远场一阶广义状态系统的有限模态。

返回 `lambda=i*kz`、`kz`、原始变量 `q`、法向导数 `qz` 以及完整状态
向量。无限特征值被保留在诊断字段 `infinite_count` 中但不作为有限
波模态返回。有限模态必须满足 `qz≈lambda*q`，并与二次远场 pencil
得到的 `kz` 一致。
"""
function solve_farfield_state_modes(base::BaseFlow;
                                    p::LSTParameters=LSTParameters(mode=:spatial),
                                    max_abs_kz::Real=100.0)
    pencil=farfield_state_pencil(base;p=p)
    E=eigen(pencil.A,pencil.E)
    vals=ComplexF64.(E.values)
    vecs=ComplexF64.(E.vectors)
    finite_ids=findall(i->isfinite(real(vals[i])) && isfinite(imag(vals[i])) &&
                              abs(-im*vals[i])<=max_abs_kz,eachindex(vals))
    sort!(finite_ids,by=i->abs(-im*vals[i]))
    modes=NamedTuple[]
    for i in finite_ids
        x=vecs[:,i]
        λ=vals[i]; q=x[1:4]; qz=x[5:8]; kz=-im*λ
        # 由广义系统直接计算残差，并单独检查导数状态关系。
        scale=max(norm(x),1.0)
        pencil_residual=norm(pencil.A*x-λ*pencil.E*x)/scale
        derivative_residual=norm(qz-λ*q)/max(norm(q),1.0)
        quadratic_residual=norm(farfield_matrix(base;p=p,kz=kz)*q)/max(norm(q),1.0)
        push!(modes,(lambda=λ,kz=kz,state=x,q=q,qz=qz,
                     pencil_residual=pencil_residual,
                     derivative_residual=derivative_residual,
                     quadratic_residual=quadratic_residual))
    end
    infinite_count=count(i->!isfinite(real(vals[i])) || !isfinite(imag(vals[i])),eachindex(vals))
    return (modes=modes,infinite_count=infinite_count,pencil=pencil,
            all_eigenvalues=vals)
end

"""
提取向 `z=+∞` 有界的有限模态子空间，并构造其左零空间约束。

在当前约定下有界条件为 `imag(kz)>imag_tolerance`。返回的 `C` 满足
`C*V≈0`，其中 `V` 是选中模态的完整状态矩阵；因此散射状态可以写成
`C*(X-X_inc)=0`。这里的 `C` 只是模态边界算子，尚未自动接入有限域
求解器，也不把有界模态等同于物理出射模态。
"""
function farfield_bounded_subspace(base::BaseFlow;
                                   p::LSTParameters=LSTParameters(mode=:spatial),
                                   imag_tolerance::Real=1e-9,
                                   max_abs_kz::Real=100.0)
    result=solve_farfield_state_modes(base;p=p,max_abs_kz=max_abs_kz)
    selected=[m for m in result.modes if imag(m.kz)>imag_tolerance]
    length(selected)>0 || error("没有找到向 z=+∞ 有界的有限远场模态")
    V=hcat([m.state for m in selected]...)
    F=svd(V;full=true)
    r=length(selected)
    size(F.U,2)>=r || error("状态模态矩阵的 SVD 维度不足")
    C=adjoint(F.U[:,r+1:end])
    return (modes=selected,V=V,C=C,rank=r,
            constraint_residual=norm(C*V),condition=cond(V),
            all_modes=result.modes,infinite_count=result.infinite_count)
end

"""
把状态约束 `C*[q;q_z]=0` 映射成有限 Chebyshev 网格上的矩阵行。

`D` 是当前有限区间的 Chebyshev 一阶微分矩阵，`edge_index` 默认为
外缘节点。返回矩阵可直接替换有限域离散系统中的外缘边界行。
"""
function farfield_state_boundary_matrix(D::AbstractMatrix,C::AbstractMatrix;
                                        edge_index::Int=size(D,1))
    n=size(D,1); size(D,2)==n || error("D 必须为方阵")
    1<=edge_index<=n || error("edge_index 超出网格范围")
    size(C,2)==8 || error("C 必须作用于 [q;q_z] 的 8 个状态分量")
    B=zeros(ComplexF64,size(C,1),4n)
    for row in 1:size(C,1), component in 1:4
        cols=(component-1)*n+1:component*n
        B[row,cols] .+= C[row,component+4].*D[edge_index,:]
        B[row,(component-1)*n+edge_index] += C[row,component]
    end
    return B
end

"""
提取奇异 descriptor pencil 的代数约束。

若 `L` 是 `E'` 的零空间基，则所有有限状态满足
`L' * A * X = 0`。这些约束属于微分代数系统内部的压力/连续性约束，
不应再次作为独立的外缘 Dirichlet 条件强行施加。
"""
function farfield_algebraic_constraints(base::BaseFlow;
                                        p::LSTParameters=LSTParameters(mode=:spatial))
    pencil=farfield_state_pencil(base;p=p)
    L=nullspace(adjoint(pencil.E))
    C=adjoint(L)*pencil.A
    finite=solve_farfield_state_modes(base;p=p)
    V=hcat([m.state for m in finite.modes]...)
    return (E=pencil.E,A=pencil.A,left_null_basis=L,C=C,
            rank=rank(C),constraint_residual=norm(C*V),
            finite_modes=finite.modes)
end

"""
诊断远场法向波数模态的压力、散度和涡量性质。

`wave_divergence` 使用局部 Cartesian 形式
`im*(alpha*u+beta*v+kz*w)`；`operator_continuity_residual` 使用当前
原始圆柱算子保留的 `Ro/R*u` 几何项。二者不同正是圆柱局部模型和
均匀 Cartesian 波模型的差别。
"""
function farfield_mode_diagnostics(base::BaseFlow;
                                   p::LSTParameters=LSTParameters(mode=:spatial),
                                   max_abs_kz::Real=100.0,
                                   kscale::Real=sqrt(p.alpha^2+p.beta^2),
                                   potential_tolerance::Real=1e-6)
    kscale>0 || error("kscale 必须为正")
    modes=solve_farfield_modes(base;p=p,max_abs_kz=max_abs_kz)
    result=NamedTuple[]
    for m in modes
        raw=ComplexF64.(m.amplitude)
        vel=raw[1:3]; velocity_norm=norm(vel)
        velocity_norm>0 || continue
        a=vel/velocity_norm; pressure=raw[4]/velocity_norm
        k=ComplexF64[p.alpha,p.beta,m.kz]
        wave_vorticity=im*cross(k,a)
        wave_divergence=im*dot(k,a)
        operator_continuity=(im*p.alpha+p.Ro/p.R)*a[1]+
                            im*p.beta*a[2]+im*m.kz*a[3]
        vorticity_ratio=norm(wave_vorticity)/kscale
        classification=vorticity_ratio<=potential_tolerance ?
                       :potential_like : :vortical_like
        push!(result,(kz=m.kz,amplitude=a,pressure=pressure,
                      pressure_ratio=abs(pressure),
                      vorticity_ratio=vorticity_ratio,
                      wave_divergence=abs(wave_divergence),
                      operator_continuity_residual=abs(operator_continuity),
                      residual=norm(m.residual)/velocity_norm,
                      classification=classification))
    end
    return result
end

"""半无限映射 z∈[0,+∞)，最后一个点表示远场。"""
function semi_infinite_grid(N::Int)
    x, Dx = _chebyshev_matrix(N)
    z=fill(Inf,N+1); scale=zeros(N+1)
    for j in eachindex(x)
        xi=x[j]
        if xi < 1-10eps(Float64)
            den=-2xi^3+xi^2-3xi+4
            z[j]=max((4xi^3-2xi^2+6xi+12)/den,0.0)
            scale[j]=(2xi^3-xi^2+3xi-4)^2/(20*(6xi^2-2xi+3))
        end
    end
    return z, Diagonal(scale)*Dx
end

_diag(x)=Diagonal(vec(x))

"""
空间二次特征值问题的块算子。

这些块对应 notebook 中的 `Ta,A,B,C,dC,D1,Vxx,Vyy,Vzz` 等算子，
并采用同一组原始变量 `(u,v,w,p)`、同一半无限 Chebyshev 网格和同一
边界约束。`L0,L1,L2` 由这些块在 `assemble_spatial` 中组装。
"""
struct SpatialOperatorBlocks
    Ta::Matrix{ComplexF64}
    A::Matrix{ComplexF64}
    B::Matrix{ComplexF64}
    C::Matrix{ComplexF64}
    dC::Matrix{ComplexF64}
    D1::Matrix{ComplexF64}
    Vxx::Matrix{ComplexF64}
    Vyy::Matrix{ComplexF64}
    Vzz::Matrix{ComplexF64}
    dVzz::Matrix{ComplexF64}
    d2Vzz::Matrix{ComplexF64}
    Vxy::Matrix{ComplexF64}
    Vxz::Matrix{ComplexF64}
    dVxz::Matrix{ComplexF64}
    Vyz::Matrix{ComplexF64}
    dVyz::Matrix{ComplexF64}
    z::Vector{Float64}
    D::Matrix{ComplexF64}
    D2::Matrix{ComplexF64}
    keep::Vector{Int}
    params::LSTParameters
end

function _sample_base(base, z)
    n=length(z); U=zeros(n); V=zeros(n); W=zeros(n)
    dU=zeros(n); dV=zeros(n); dW=zeros(n)
    for j in 1:n
        U[j],V[j],W[j],dU[j],dV[j],dW[j]=_profile_at(base,z[j])
    end
    return U,V,W,dU,dV,dW
end

"""组装时间模式 A*q=omega*B*q。"""
function assemble_temporal(base::BaseFlow; p::LSTParameters=LSTParameters())
    N=p.N; n=N+1; z,D=finite_grid(N,p.temporal_zmax); D2=D*D
    U,V,W,dU,dV,dW=_sample_base(base,z); ddU=D*dU; ddV=D*dV
    Id=Matrix{Float64}(I,n,n); α=p.alpha; β=p.beta; R=p.R; Ro=p.Ro; Co=p.Co
    αb=α-im/R; λ2=α^2+β^2; λb2=α*αb+β^2
    L=D2-λ2*Id; Lb=D2-λb2*Id
    Ud=_diag(U); Vd=_diag(V); Wd=_diag(W); dUd=_diag(dU); dVd=_diag(dV)
    dWd=_diag(dW); ddUd=_diag(ddU); ddVd=_diag(ddV); Cd=_diag(α.*U.+β.*V)
    A11=im*(L*Lb)+R*Cd*Lb-R*(αb*ddUd+β*ddVd)-im*Ro*Wd*D*Lb-im*Ro*dWd*Lb-im*Ro*Ud*D2
    A12=_diag(2Ro.*V.+Co)*D+2Ro*dVd
    A21=_diag(2Ro.*V.+Co)*D-im*R*_diag(α.*dV.-β.*dU)
    A22=im*L+R*Cd-im*Ro*Wd*D-im*Ro*Ud
    B11=R*Lb; B22=R*Matrix{ComplexF64}(I,n,n)
    A=ComplexF64.([A11 A12; A21 A22])
    B=ComplexF64.([B11 zeros(ComplexF64,n,n); zeros(ComplexF64,n,n) B22])
    rows=[1,2,n-1,n,n+1,2n]
    for r in rows; A[r,:].=0; B[r,:].=0; end
    A[1,1]=1; A[2,1:n].=D[1,:]; A[n-1,n]=1; A[n,1:n].=D[n,:]
    A[n+1,n+1]=1; A[2n,2n]=1
    return (A=A,B=B,z=z,D=D,params=p)
end

function solve_temporal(op; max_abs_omega=1.0)
    E=eigen(op.A,op.B); vals=ComplexF64.(E.values); vecs=ComplexF64.(E.vectors)
    ids=findall(i->isfinite(real(vals[i]))&&isfinite(imag(vals[i]))&&abs(vals[i])<=max_abs_omega,eachindex(vals))
    sort!(ids,by=i->abs(vals[i]))
    return [(omega=vals[i],vector=vecs[:,i]) for i in ids]
end

function _block4(a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44)
    return [a11 a12 a13 a14; a21 a22 a23 a24; a31 a32 a33 a34; a41 a42 a43 a44]
end

"""
组装 notebook 风格的空间块算子。

块算子只保存原始边界条件施加前的矩阵。空间二次铅笔为
`L0 + alpha*L1 + alpha^2*L2 = 0`，其中 `L0,L1,L2` 在
`assemble_spatial` 中从这些块组装并删去同一组边界自由度。
`grid=:affine` 用于具有非零有限外缘输入的受迫问题。
"""
function assemble_spatial_blocks(base::BaseFlow;
                                 p::LSTParameters=LSTParameters(mode=:spatial),
                                 grid::Symbol=:semi_infinite,
                                 zmax::Float64=30.0)
    N=p.N; n=N+1
    z,D0 = if grid == :semi_infinite
        semi_infinite_grid(N)
    elseif grid == :finite
        finite_grid(N,zmax)
    elseif grid == :affine
        affine_grid(N,zmax)
    else
        error("grid 必须是 :semi_infinite、:finite 或 :affine")
    end
    D20=D0*D0
    F,G,H,_,_,_=_sample_base(base,z); R=p.R; Ro=p.Ro; Co=p.Co
    Id=Matrix{ComplexF64}(I,n,n); Z=zeros(ComplexF64,n,n)
    D=Matrix{ComplexF64}(D0); D2=Matrix{ComplexF64}(D20)
    Fd=Diagonal(ComplexF64.(F)); Gd=Diagonal(ComplexF64.(G)); Hd=Diagonal(ComplexF64.(H))
    # 原始 BEK 算子中的 D*F、D*G、D*H 是系数的法向导数，写成对角矩阵。
    Dfd=Diagonal(ComplexF64.(D*F)); Dgd=Diagonal(ComplexF64.(D*G))
    Dhd=(Ro/R)*Diagonal(ComplexF64.(D*H))
    Cor=Diagonal(ComplexF64.(2Ro.*G.+Co))

    Ta=[Id Z Z Z; Z Id Z Z; Z Z Id Z; Z Z Z Z]
    A=[Fd Z Z Id; Z Fd Z Z; Z Z Fd Z; Id Z Z Z]
    B=[(1/R)*Gd Z Z Z; Z (1/R)*Gd Z (1/R)*Id; Z Z (1/R)*Gd Z; Z (1/R)*Id Z Z]
    C=[(Ro/R)*Hd Z Z Z; Z (Ro/R)*Hd Z Z; Z Z (Ro/R)*Hd Id; Z Z Id Z]
    dC=[Dhd Z Z Z; Z Dhd Z Z; Z Z Dhd Z; Z Z Z Z]
    D1=[(Ro/R)*Fd -(1/R)*Cor Dfd Z;
        (1/R)*Cor (Ro/R)*Fd Dgd Z;
        Z Z Dhd Z;
        (Ro/R)*Id Z Z Z]
    Vxx=[-(1/R)*Id Z Z Z; Z -(1/R)*Id Z Z; Z Z -(1/R)*Id Z; Z Z Z Z]
    Vyy=[-(1/R^3)*Id Z Z Z; Z -(1/R^3)*Id Z Z; Z Z -(1/R^3)*Id Z; Z Z Z Z]
    Vzz=[-(1/R)*Id Z Z Z; Z -(1/R)*Id Z Z; Z Z -(1/R)*Id Z; Z Z Z Z]
    # 当前不可压 von Kármán 算子没有横向黏性耦合块；保留这些字段以便
    # 与 receptivity.ipynb 的 A/B/C/D 算子接口逐项对应。
    dVzz=zeros(ComplexF64,4n,4n); d2Vzz=zeros(ComplexF64,4n,4n)
    Vxy=zeros(ComplexF64,4n,4n); Vxz=zeros(ComplexF64,4n,4n)
    dVxz=zeros(ComplexF64,4n,4n); Vyz=zeros(ComplexF64,4n,4n)
    dVyz=zeros(ComplexF64,4n,4n)
    bc=[1,n,n+1,2n,2n+1,3n,4n]
    keep=collect(setdiff(1:4n,bc))
    return SpatialOperatorBlocks(Ta,A,B,C,dC,D1,Vxx,Vyy,Vzz,dVzz,d2Vzz,
                                Vxy,Vxz,dVxz,Vyz,dVyz,z,D,D2,keep,p)
end

"""
组装二次原始变量空间算子 L0+alpha*L1+alpha^2*L2=0。

在 Malik 的等温 Type-I 参数下，这个原始变量离散问题的收敛值约为
0.38505249+3.48666e-4im；它与消压后的文献中性点
0.38482+0im 存在系统差异。因此这里的残差和网格收敛不能表述为严格复现
Malik 的 alpha，后者还取决于远场消压渐近条件和局部方程截断约定。
"""
function assemble_spatial(base::BaseFlow; p::LSTParameters=LSTParameters(mode=:spatial))
    blocks=assemble_spatial_blocks(base;p=p); n=p.N+1; keep=blocks.keep
    L0,L1,L2=_spatial_pencil(blocks,p)
    return (L0=Matrix{ComplexF64}(L0[keep,keep]),
            L1=Matrix{ComplexF64}(L1[keep,keep]),
            L2=Matrix{ComplexF64}(L2[keep,keep]),
            z=blocks.z,D=blocks.D,params=p,blocks=blocks)
end

"""从块算子生成未施加边界条件的空间二次铅笔。"""
function _spatial_pencil(blocks::SpatialOperatorBlocks,p::LSTParameters)
    Kd=kron(Matrix{ComplexF64}(I,4,4),blocks.D)
    Kd2=kron(Matrix{ComplexF64}(I,4,4),blocks.D2)
    L0=blocks.D1+im*p.beta*p.R*blocks.B-im*p.omega*blocks.Ta-
       p.beta^2*p.R^2*blocks.Vyy+blocks.C*Kd+blocks.Vzz*Kd2
    L1=im*blocks.A
    L2=-blocks.Vxx
    return L0,L1,L2
end

"""有限外缘上的提升函数形状。"""
function _lift_shape(s::Float64,lift::Symbol)
    ξ=clamp(s,0.0,1.0)
    if lift == :linear
        return ξ
    elseif lift == :smooth
        return ξ^2*(3.0-2.0*ξ)
    else
        error("lift 必须是 :linear 或 :smooth")
    end
end

"""构造满足壁面零值、外缘给定值的四变量计算提升场。"""
function _forced_lift(z::Vector{Float64},zmax::Float64,
                      edge_state,lift::Symbol)
    length(edge_state)==4 || error("edge_state 必须包含 [u_e,v_e,w_e,p_e]")
    n=length(z); edge=ComplexF64.(collect(edge_state)); ell=zeros(ComplexF64,4n)
    for j in 1:n
        g=_lift_shape(z[j]/zmax,lift)
        ell[j]=g*edge[1]
        ell[n+j]=g*edge[2]
        ell[2n+j]=g*edge[3]
        ell[3n+j]=g*edge[4]
    end
    return ell
end

"""组装固定外部 alpha 的有限外缘受迫空间矩阵。"""
function assemble_forced_spatial(base::BaseFlow;
                                 p::LSTParameters=LSTParameters(mode=:spatial),
                                 zmax::Float64=30.0)
    zmax>0 || error("zmax 必须为正")
    blocks=assemble_spatial_blocks(base;p=p,grid=:affine,zmax=zmax)
    L0,L1,L2=_spatial_pencil(blocks,p)
    P=L0+p.alpha*L1+p.alpha^2*L2
    return (P=Matrix{ComplexF64}(P),L0=L0,L1=L1,L2=L2,
            z=blocks.z,D=blocks.D,params=p,blocks=blocks,zmax=zmax)
end

"""
求解光滑盘面对给定外缘涡波的直接线性响应。

`edge_state=[u_e,v_e,w_e,p_e]` 是有限外缘 z=zmax 处的复幅值。
提升场只用于把非齐次外缘条件转成齐次修正问题，最终返回的 `q` 是
恢复后的总扰动。`source` 若给出，应是完整 4*(N+1) 向量形式的额外
体源；默认没有粗糙元，因此 source 默认为零。`wall_velocity` 是壁面
处的复速度边界 `[u_w,v_w,w_w]`，默认是无滑移的零向量。压力不在壁面
指定；外缘压力仍由 `edge_state[4]` 给出。
"""
function solve_forced_spatial(base::BaseFlow,edge_state;
                              p::LSTParameters=LSTParameters(mode=:spatial),
                              zmax::Float64=30.0,
                              lift::Symbol=:smooth,
                              wall_velocity=zeros(ComplexF64,3),
                              singular_strategy::Symbol=:error,
                              source=nothing)
    op=assemble_forced_spatial(base;p=p,zmax=zmax)
    blocks=op.blocks; n=p.N+1; P=op.P
    ell=_forced_lift(op.z,zmax,edge_state,lift)
    src=source===nothing ? zeros(ComplexF64,4n) : ComplexF64.(collect(source))
    length(src)==4n || error("source 必须有 4*(N+1) 个分量")
    wall=ComplexF64.(collect(wall_velocity))
    length(wall)==3 || error("wall_velocity 必须包含 [u_w,v_w,w_w]")
    # 提升函数产生的等效体源；这不是额外的物理体力。
    rhs=src-P*ell
    bc=collect(setdiff(1:4n,blocks.keep))
    boundary_values=zeros(ComplexF64,4n)
    edge=ComplexF64.(collect(edge_state))
    boundary_values[1]=wall[1]
    boundary_values[n]=edge[1]
    boundary_values[n+1]=wall[2]
    boundary_values[2n]=edge[2]
    boundary_values[2n+1]=wall[3]
    boundary_values[3n]=edge[3]
    boundary_values[4n]=edge[4]
    Pbc=copy(P); rhs_bc=copy(rhs)
    # 对修正场施加总场的目标边界值；压力只在外缘指定。
    for row in bc
        Pbc[row,:].=0
        Pbc[row,row]=1
        rhs_bc[row]=boundary_values[row]-ell[row]
    end
    correction=if singular_strategy==:error
        Pbc\rhs_bc
    elseif singular_strategy==:pinv
        # 仅用于 κ=0 等已诊断的压力零模；速度场需另行检查。
        pinv(Pbc;rtol=1e-12)*rhs_bc
    else
        error("singular_strategy 必须是 :error 或 :pinv")
    end
    q=ell+correction
    residual=P*q-src
    interior_residual=residual[blocks.keep]
    boundary_error=q[bc]-boundary_values[bc]
    qz=vcat(op.D*q[1:n],op.D*q[n+1:2n],op.D*q[2n+1:3n],op.D*q[3n+1:4n])
    return (q=q,correction=correction,lift=ell,source=src,
            edge_state=edge, P=P,P_boundary=Pbc,rhs=rhs_bc,
            residual=residual,interior_residual=interior_residual,
            boundary_error=boundary_error,z=op.z,D=op.D,
            qz=qz,velocity_derivatives=qz[1:3n],
            params=p,blocks=blocks,zmax=zmax,lift_kind=lift)
end

"""从内部 von Kármán ODE 基流直接求解光滑盘面受迫响应。"""
function solve_forced_spatial(edge_state;
                              p::LSTParameters=LSTParameters(mode=:spatial),
                              base_params::VonKarmanBaseParameters=VonKarmanBaseParameters(),
                              zmax::Float64=30.0,
                              lift::Symbol=:smooth,
                              wall_velocity=zeros(ComplexF64,3),
                              singular_strategy::Symbol=:error,
                              source=nothing)
    base=solve_von_karman_baseflow(params=base_params)
    return solve_forced_spatial(base,edge_state;p=p,zmax=zmax,lift=lift,
                                wall_velocity=wall_velocity,
                                singular_strategy=singular_strategy,source=source)
end

"""计算光滑盘响应在壁面的复速度梯度。

输入可以是 `solve_forced_spatial` 返回的响应对象，也可以显式给出
`q` 和 `D`。节点的第一个点是 z=0，因此返回
`[u_z(0),v_z(0),w_z(0)]`。压力梯度不属于壁面速度散射边界条件。
"""
function wall_shear_response(q::AbstractVector,D::AbstractMatrix)
    n=size(D,1); n>=2 || error("响应网格至少需要两个节点")
    length(q)==4n || error("q 必须包含 4*(N+1) 个分量")
    size(D)==(n,n) || error("D 必须是与响应网格一致的微分矩阵")
    gradients=ComplexF64[]
    for c in 1:3
        values=ComplexF64.(q[(c-1)*n+1:c*n])
        push!(gradients,(D*values)[1])
    end
    return gradients
end

function wall_shear_response(response; q=response.q, D=response.D)
    if hasproperty(response,:wall_shear)
        return ComplexF64.(collect(getproperty(response,:wall_shear)))
    end
    return wall_shear_response(q,D)
end

"""在任意壁面速度 Dirichlet 条件下求解有限外缘空间响应。

这是 `solve_forced_spatial` 的语义化接口，`wall_velocity` 直接对应
z=0 处的 `[u,v,w]`，`edge_state` 仍对应 z=zmax 处的 `[u,v,w,p]`。
"""
function solve_wall_forced_spatial(base::BaseFlow,edge_state;
                                   wall_velocity,
                                   p::LSTParameters=LSTParameters(mode=:spatial),
                                   zmax::Float64=30.0,
                                   lift::Symbol=:smooth,
                                   singular_strategy::Symbol=:error,
                                   source=nothing)
    return solve_forced_spatial(base,edge_state;p=p,zmax=zmax,lift=lift,
                                wall_velocity=wall_velocity,
                                singular_strategy=singular_strategy,source=source)
end

"""使用内部 von Kármán 基流求解器计算任意壁面速度响应。"""
function solve_wall_forced_spatial(edge_state;
                                   wall_velocity,
                                   p::LSTParameters=LSTParameters(mode=:spatial),
                                   base_params::VonKarmanBaseParameters=VonKarmanBaseParameters(),
                                   zmax::Float64=30.0,
                                   lift::Symbol=:smooth,
                                   singular_strategy::Symbol=:error,
                                   source=nothing)
    base=solve_von_karman_baseflow(params=base_params)
    return solve_wall_forced_spatial(base,edge_state;wall_velocity=wall_velocity,
                                     p=p,zmax=zmax,lift=lift,
                                     singular_strategy=singular_strategy,source=source)
end

"""构造 Gaussian 径向包络及其周向 Fourier 分支。

`branch=:real` 返回 `epsilon*H_r(r)*cos(n*theta)`；`branch=:plus`
和 `branch=:minus` 分别返回复分支 `epsilon*H_r(r)*exp(±im*n*theta)`。
这里的 `epsilon` 是物理粗糙度幅值，单位 Fourier 分支本身应取
`epsilon=1` 后再在重构时统一乘入。
"""
function gaussian_roughness(r;
                           rf::Real,
                           c::Real,
                           theta::Real=0.0,
                           n::Real=0.0,
                           epsilon::Number=1.0,
                           branch::Symbol=:real)
    c>0 || error("Gaussian 宽度 c 必须为正")
    Hr=exp.(-((r.-rf).^2)./(2c^2))
    amp=ComplexF64(epsilon)
    if branch==:real
        return real.(amp.*Hr.*cos(n*theta))
    elseif branch==:plus
        return amp.*Hr.*exp(im*n*theta)
    elseif branch==:minus
        return amp.*Hr.*exp(-im*n*theta)
    else
        error("branch 必须是 :real、:plus 或 :minus")
    end
end

"""Gaussian 径向包络的解析 Fourier 谱。

采用

    Hhat(kappa)=∫ H_r(r)exp(-im*kappa*r)dr,
    H_r(r)=(1/(2π))*∫ Hhat(kappa)exp(im*kappa*r)dkappa.

因此 `Hhat=sqrt(2π)*c*exp(-c^2*kappa^2/2)*exp(-im*kappa*rf)`。
"""
function roughness_fourier_spectrum(kappa;
                                    rf::Real,
                                    c::Real)
    c>0 || error("Gaussian 宽度 c 必须为正")
    return sqrt(2π)*c .* exp.(-0.5*c^2 .* kappa.^2) .* exp.(-im .* kappa .* rf)
end

"""生成 Gaussian 频谱的均匀 κ 网格和梯形积分权重。"""
function roughness_fourier_grid(; c::Real,
                                kappa_max::Real=6.0/c,
                                nk::Int=128)
    c>0 || error("Gaussian 宽度 c 必须为正")
    kappa_max>0 || error("kappa_max 必须为正")
    nk>=3 || error("nk 至少取 3")
    kappa=collect(range(-Float64(kappa_max),Float64(kappa_max);length=nk))
    dk=kappa[2]-kappa[1]
    weights=fill(dk,nk); weights[[1,end]].*=0.5
    return (kappa=kappa,weights=weights,kappa_max=Float64(kappa_max),nk=nk)
end

function _state_vector(q)
    if q isa AbstractVector
        return ComplexF64.(collect(q))
    elseif hasproperty(q,:q)
        return ComplexF64.(collect(getproperty(q,:q)))
    else
        error("谱响应必须是 4n 向量或包含 q 字段的响应对象")
    end
end

function _reconstruct_roughness_branch(kappa,qhat;
                                      r::Real,
                                      theta::Real,
                                      n::Real,
                                      epsilon::Number,
                                      rf::Real,
                                      c::Real,
                                      weights)
    length(kappa)==length(qhat) || error("kappa 与 qhat 长度必须一致")
    length(kappa)>=2 || error("Fourier 重构至少需要两个 κ 点")
    kk=Float64.(collect(kappa)); dk=diff(kk)
    any(dk.<=0) && error("kappa 必须严格递增")
    if weights===nothing
        d0=sum(dk)/length(dk)
        maximum(abs.(dk.-d0))<=1e-10*max(1.0,abs(d0)) ||
            error("非均匀 kappa 网格必须显式提供积分 weights")
        w=fill(d0,length(kk)); w[[1,end]].*=0.5
    else
        w=Float64.(collect(weights))
    end
    length(w)==length(kappa) || error("weights 与 kappa 长度必须一致")
    a0=_state_vector(qhat[1]); out=zeros(ComplexF64,length(a0))
    for j in eachindex(kappa)
        aj=_state_vector(qhat[j]); length(aj)==length(a0) ||
            error("所有 qhat 必须具有相同长度")
        coef=ComplexF64(epsilon)*w[j]/(2π)*
             roughness_fourier_spectrum(kk[j];rf=rf,c=c)*exp(im*kk[j]*r+im*n*theta)
        out .+= coef.*aj
    end
    return out
end

"""根据单位 Fourier 响应重构 Gaussian 粗糙度场。

当 `branch=:plus` 时，`qhat` 是 (+n) 分支；当 `branch=:minus` 时，
`qhat` 是 (-n) 分支。`branch=:real` 时，`qhat` 必须是包含 `plus`
和 `minus` 两个字段的 NamedTuple；对 n≠0 自动给两个分支各乘 1/2，
而 n=0 只保留一个分支。`qhat[j]` 应是单位粗糙度幅值下第 j 个 κ 的
响应，`epsilon` 只在本函数中乘一次。
"""
function reconstruct_roughness_field(kappa,qhat;
                                    r::Real,
                                    theta::Real=0.0,
                                    n::Real=0.0,
                                    epsilon::Number=1.0,
                                    rf::Real,
                                    c::Real,
                                    weights=nothing,
                                    branch::Symbol=:plus)
    if branch==:real
        hasproperty(qhat,:plus) && hasproperty(qhat,:minus) ||
            error("branch=:real 时 qhat 必须含 plus 和 minus 字段")
        plus=_reconstruct_roughness_branch(kappa,qhat.plus;r=r,theta=theta,
                                           n=Float64(n),epsilon=epsilon,rf=rf,c=c,
                                           weights=weights)
        if iszero(n)
            return plus
        end
        minus=_reconstruct_roughness_branch(kappa,qhat.minus;r=r,theta=theta,
                                            n=-Float64(n),epsilon=epsilon,rf=rf,c=c,
                                            weights=weights)
        return 0.5 .* (plus.+minus)
    elseif branch==:plus
        return _reconstruct_roughness_branch(kappa,qhat;r=r,theta=theta,n=Float64(n),
                                             epsilon=epsilon,rf=rf,c=c,weights=weights)
    elseif branch==:minus
        return _reconstruct_roughness_branch(kappa,qhat;r=r,theta=theta,n=-Float64(n),
                                             epsilon=epsilon,rf=rf,c=c,weights=weights)
    else
        error("branch 必须是 :plus、:minus 或 :real")
    end
end

"""重构交叉阶响应；接口与 `reconstruct_roughness_field` 相同。"""
function reconstruct_coupled_field(kappa,qhat; kwargs...)
    return reconstruct_roughness_field(kappa,qhat;kwargs...)
end

"""单一 Fourier 粗糙度分量产生的壁面散射响应。

`smooth_response` 是光滑盘面对外部涡波的响应；粗糙度高度为
`h_amplitude*exp(i*(kappa*r+n*theta))`。局部 Fourier 测试中，目标
波数取 `alpha_target=alpha_f+kappa`；周向目标波数必须通过
`target_beta` 或 `beta_shift` 显式给出，`n`
只作为 Fourier 分支标签。函数只构造壁面 Taylor 项，
不包含粗糙基本流畸变导致的体积交互强迫。
"""
function roughness_wall_scattering(base::BaseFlow,smooth_response;
                                   edge_state,
                                   p::LSTParameters=LSTParameters(mode=:spatial),
                                   kappa::Real=0.0,
                                   n::Real=0.0,
                                   h_amplitude::Number=1e-3,
                                   zmax::Float64=30.0,
                                   lift::Symbol=:smooth,
                                   singular_strategy::Symbol=:error,
                                   beta_shift::Union{Nothing,Real}=nothing,
                                   target_alpha::Union{Nothing,Real}=nothing,
                                   target_beta::Union{Nothing,Real}=nothing,
                                   source=nothing)
    shear=wall_shear_response(smooth_response)
    wall_velocity=-ComplexF64(h_amplitude).*ComplexF64.(shear)
    alpha_target=target_alpha===nothing ? p.alpha+Float64(kappa) : Float64(target_alpha)
    beta_target=if target_beta!==nothing
        Float64(target_beta)
    elseif beta_shift!==nothing
        p.beta+Float64(beta_shift)
    elseif iszero(n)
        p.beta
    else
        error("n 非零时必须显式提供 target_beta 或 beta_shift；当前代码没有确认整数 m 到局部 beta 的换算")
    end
    pt=LSTParameters(mode=:spatial,R=p.R,alpha=alpha_target,beta=beta_target,
                     omega=p.omega,Ro=p.Ro,Co=p.Co,N=p.N,
                     temporal_zmax=p.temporal_zmax,spatial_shift=p.spatial_shift,
                     nev=p.nev)
    response=solve_wall_forced_spatial(base,edge_state;wall_velocity=wall_velocity,
                                       p=pt,zmax=zmax,lift=lift,
                                       singular_strategy=singular_strategy,source=source)
    return (response=response,wall_velocity=wall_velocity,shear=shear,
            h_amplitude=ComplexF64(h_amplitude),kappa=Float64(kappa),
            roughness_mode=Float64(n),params=pt)
end

"""求解单一 Fourier 粗糙度引起的一阶定常基本流畸变。

粗糙度取局部形式 `h_amplitude*exp(i*(kappa*r+n*theta))`，并令
`omega=0`。壁面 Taylor 条件为

    U_h(0) = -h_amplitude*[U_0'(0), V_0'(0), W_0'(0)]

默认外缘速度和压力均为零，这对应有限计算区间内的局部粗糙度原型。
`n` 在当前局部算子中按物理周向波数使用；若输入的是整数周向模态，
应在调用前完成半径换算。该函数只求解几何壁面强迫，不包含后续
涡波–粗糙度交叉阶体积源项。
"""
function solve_roughness_distortion(base::BaseFlow;
                                    p::LSTParameters=LSTParameters(mode=:spatial),
                                    kappa::Real=0.1,
                                    n::Real=0.0,
                                    h_amplitude::Number=1.0,
                                    zmax::Float64=30.0,
                                    lift::Symbol=:smooth,
                                    singular_strategy::Symbol=:error,
                                    edge_state=zeros(ComplexF64,4),
                                    beta_shift::Union{Nothing,Real}=nothing,
                                    target_alpha::Union{Nothing,Real}=nothing,
                                    target_beta::Union{Nothing,Real}=nothing,
                                    source=nothing)
    alpha_h=target_alpha===nothing ? Float64(kappa) : Float64(target_alpha)
    beta_h=if target_beta!==nothing
        Float64(target_beta)
    elseif beta_shift!==nothing
        Float64(beta_shift)
    elseif iszero(n)
        0.0
    else
        error("n 非零时必须显式提供 target_beta 或 beta_shift；当前代码没有确认整数 m 到局部 beta 的换算")
    end
    ph=LSTParameters(mode=:spatial,R=p.R,alpha=alpha_h,beta=beta_h,omega=0.0,
                     Ro=p.Ro,Co=p.Co,N=p.N,temporal_zmax=p.temporal_zmax,
                     spatial_shift=p.spatial_shift,nev=p.nev)
    if iszero(alpha_h) && iszero(beta_h) && singular_strategy==:error
        error("κ=0、beta_h=0 的定常均匀 Fourier 分支含离散压力零模；请显式使用 singular_strategy=:pinv 并检查速度场，或避开该点")
    end
    baseflow_shear=ComplexF64[base.du[1],base.dv[1],base.dw[1]]
    wall_velocity=-ComplexF64(h_amplitude).*baseflow_shear
    response=solve_wall_forced_spatial(base,edge_state;
                                       wall_velocity=wall_velocity,p=ph,
                                       zmax=zmax,lift=lift,
                                       singular_strategy=singular_strategy,source=source)
    return (response=response,wall_velocity=wall_velocity,
            baseflow_shear=baseflow_shear,h_amplitude=ComplexF64(h_amplitude),
            kappa=Float64(kappa),roughness_mode=Float64(n),params=ph)
end

"""使用内部 von Kármán 基流求解器计算一阶定常粗糙度畸变。"""
function solve_roughness_distortion(; p::LSTParameters=LSTParameters(mode=:spatial),
                                    base_params::VonKarmanBaseParameters=VonKarmanBaseParameters(),
                                    kwargs...)
    base=solve_von_karman_baseflow(params=base_params)
    return solve_roughness_distortion(base;p=p,kwargs...)
end

function _state_data(q; D=nothing)
    qv=_state_vector(q); length(qv)%4==0 || error("状态向量长度必须是 4n")
    n=length(qv)÷4
    Dm=if D!==nothing
        Matrix{ComplexF64}(D)
    elseif hasproperty(q,:D)
        Matrix{ComplexF64}(getproperty(q,:D))
    else
        error("输入状态不是响应对象时必须显式提供 D")
    end
    size(Dm)==(n,n) || error("D 的尺寸必须与状态向量一致")
    qzv=if hasproperty(q,:qz)
        z=ComplexF64.(collect(getproperty(q,:qz)))
        length(z)==4n || error("qz 必须包含 4n 个分量")
        z
    else
        nothing
    end
    return (q=qv,n=n,D=Dm,qz=qzv)
end

function _local_bilinear_convection(a,b,pa,pb,D,curv,n; b_qz=nothing)
    ua=a[1:n]; va=a[n+1:2n]; wa=a[2n+1:3n]
    ub=b[1:n]; vb=b[n+1:2n]; wb=b[2n+1:3n]
    adv=(im*pb.alpha).*ua .+ (im*pb.beta).*va
    dbu=b_qz===nothing ? D*ub : b_qz[1:n]
    dbv=b_qz===nothing ? D*vb : b_qz[n+1:2n]
    dbw=b_qz===nothing ? D*wb : b_qz[2n+1:3n]
    Nu=adv.*ub .+ curv.*wa.*dbu .- curv.*va.*vb
    Nv=adv.*vb .+ curv.*wa.*dbv .+ curv.*va.*ub
    Nw=adv.*wb .+ curv.*wa.*dbw
    return Nu,Nv,Nw
end

"""构造 frozen-R 圆柱局部模型中的二阶交叉体积强迫。

输入 `qh` 和 `qf` 可以是 4n 状态向量或响应对象。其参数分别由
`p_h`、`p_f` 给出；若省略，则从响应对象的 `params` 字段读取。当前
局部算子对应的双线性对流项采用

    B_r(a,b)=a_u*∂r b_u+a_v*∂s b_u+(Ro/R)*a_w*∂z b_u
             -(Ro/R)*a_v*b_v,
    B_s(a,b)=a_u*∂r b_v+a_v*∂s b_v+(Ro/R)*a_w*∂z b_v
             +(Ro/R)*a_v*b_u,
    B_z(a,b)=a_u*∂r b_w+a_v*∂s b_w+(Ro/R)*a_w*∂z b_w.

这里 `∂r=im*alpha`、`∂s=im*beta` 沿用当前 `im*alpha*A` 和
`im*beta*R*B` 的局部 Fourier 约定；压力和连续性不产生非线性源，
Coriolis 项也不进入源项。返回的 `source` 已包含交叉项的负号，
可以直接传给 `solve_forced_spatial(source=...)`。
"""
function nonlinear_cross_forcing(qh,qf,base::BaseFlow;
                                 p_h=nothing,
                                 p_f=nothing,
                                 D=nothing)
    hd=_state_data(qh;D=D); fd=_state_data(qf;D=D)
    hd.n==fd.n || error("qh 与 qf 必须使用相同的法向离散")
    maximum(abs,hd.D-fd.D)<=1e-10 || error("qh 与 qf 的微分矩阵必须一致")
    ph=p_h===nothing ?
       (hasproperty(qh,:params) ? getproperty(qh,:params) :
        error("请提供 p_h 或使用包含 params 的 qh 响应对象")) : p_h
    pf=p_f===nothing ?
       (hasproperty(qf,:params) ? getproperty(qf,:params) :
        error("请提供 p_f 或使用包含 params 的 qf 响应对象")) : p_f
    abs(ph.R-pf.R)<=1e-10 || error("qh 与 qf 必须使用相同的 R")
    abs(ph.Ro-pf.Ro)<=1e-10 || error("qh 与 qf 必须使用相同的 Ro")
    abs(ph.Co-pf.Co)<=1e-10 || error("qh 与 qf 必须使用相同的 Co")
    n=hd.n; curv=ph.Ro/ph.R
    N1=_local_bilinear_convection(hd.q,fd.q,ph,pf,hd.D,curv,n;b_qz=fd.qz)
    N2=_local_bilinear_convection(fd.q,hd.q,pf,ph,hd.D,curv,n;b_qz=hd.qz)
    cross_u=N1[1].+N2[1]; cross_v=N1[2].+N2[2]; cross_w=N1[3].+N2[3]
    source=vcat(-cross_u,-cross_v,-cross_w,zeros(ComplexF64,n))
    return (source=source,velocity=ComplexF64.([-cross_u;-cross_v;-cross_w]),
            convective=ComplexF64.([cross_u;cross_v;cross_w]),
            target_alpha=ph.alpha+pf.alpha,
            target_beta=ph.beta+pf.beta,
            target_omega=ph.omega+pf.omega,
            params_h=ph,params_f=pf,D=hd.D)
end

"""求解同时包含二阶体积源和二阶壁面几何强迫的耦合响应。

`qh` 应是单位粗糙度 Fourier 分支的定常畸变，`qf` 应是给定外缘
Dirichlet 幅值下的完整一阶响应。`roughness_amplitude` 是该 Fourier
粗糙度分支的实际复系数；Gaussian 重构时它应包含
`epsilon*Hhat(kappa)*Δkappa/(2π)` 以及 n≠0 时的 1/2 分支系数。
函数同时返回总响应、仅体积强迫响应和仅壁面强迫响应。
"""
function solve_coupled_response(base::BaseFlow,qh,qf;
                                p_h=nothing,
                                p_f=nothing,
                                roughness_amplitude::Number=1.0+0im,
                                edge_state=zeros(ComplexF64,4),
                                zmax::Float64=30.0,
                                lift::Symbol=:smooth,
                                decompose::Bool=true)
    length(edge_state)==4 || error("交叉阶 edge_state 必须包含四个分量")
    maximum(abs,ComplexF64.(collect(edge_state)))<=1e-14 ||
        error("当前二阶耦合模型要求齐次外缘条件；真实散射外缘条件尚未实现")
    forcing=nonlinear_cross_forcing(qh,qf,base;p_h=p_h,p_f=p_f)
    ph=forcing.params_h; pf=forcing.params_f
    pt=LSTParameters(mode=:spatial,R=pf.R,
                     alpha=forcing.target_alpha,beta=forcing.target_beta,
                     omega=forcing.target_omega,Ro=pf.Ro,Co=pf.Co,
                     N=pf.N,temporal_zmax=pf.temporal_zmax,
                     spatial_shift=pf.spatial_shift,nev=pf.nev)
    scale=ComplexF64(roughness_amplitude)
    source=scale.*forcing.source
    shear=wall_shear_response(qf)
    wall_velocity=-scale.*shear
    zero_wall=zeros(ComplexF64,3); zero_source=zeros(ComplexF64,length(source))
    total=solve_wall_forced_spatial(base,edge_state;wall_velocity=wall_velocity,
                                    p=pt,zmax=zmax,lift=lift,source=source)
    volume=wall= nothing
    if decompose
        # 分解项只包含二阶物理源和壁面强迫，外缘必须取齐次值，
        # 否则同一个一阶外缘 Dirichlet 输入会在两部分中重复计入。
        zero_edge=zeros(ComplexF64,4)
        volume=solve_wall_forced_spatial(base,zero_edge;wall_velocity=zero_wall,
                                         p=pt,zmax=zmax,lift=lift,source=source)
        wall=solve_wall_forced_spatial(base,zero_edge;wall_velocity=wall_velocity,
                                       p=pt,zmax=zmax,lift=lift,source=zero_source)
    end
    return (response=total,volume_response=volume,wall_response=wall,
            source=source,unit_source=forcing.source,wall_velocity=wall_velocity,
            unit_wall_velocity=-shear,forcing=forcing,params=pt,
            roughness_amplitude=scale)
end

"""计算一维节点上的 barycentric 插值权重。"""
function _barycentric_weights(x::AbstractVector{<:Real})
    n=length(x); w=ones(Float64,n)
    for j in 1:n
        for k in 1:n
            j==k && continue
            w[j]/=(x[j]-x[k])
        end
    end
    w./=maximum(abs,w)
    return w
end

"""在任意查询点评价谱离散的复值剖面。"""
function _barycentric_interpolate(x::Vector{Float64},y::Vector{ComplexF64},
                                  bw::Vector{Float64},xq::Float64)
    j=findfirst(j->abs(xq-x[j])<=32eps(Float64)*max(1.0,abs(xq)),eachindex(x))
    j!==nothing && return y[j]
    den=0.0; num=0.0+0im
    for k in eachindex(x)
        t=bw[k]/(xq-x[k])
        den+=t; num+=t*y[k]
    end
    return num/den
end

"""
计算物理入射场扣除后的响应增益。

`response.q` 是总场，`edge_state` 定义入射平面涡波在 z=zmax 处的
复幅值。`kz` 用于入射场的法向相位；kz=0 时入射场在法向保持常值。
输出积分在重新生成的 [0,zi] Clenshaw–Curtis 节点上完成。
"""
function response_gains(response;
                        edge_state,
                        zi::Float64=10.0,
                        Nq::Int=128,
                        kz::Real=0.0)
    n=length(response.z); n>=2 || error("响应网格至少需要两个节点")
    zedge=response.z[end]
    zi<=zedge || error("zi 不能超过响应外缘 zmax")
    length(edge_state)==4 || error("edge_state 必须包含四个分量")
    edge=ComplexF64.(collect(edge_state))
    zq,wq=subinterval_cc_grid(zi,Nq)
    bw=_barycentric_weights(response.z)
    total=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        values=ComplexF64.(response.q[(c-1)*n+1:c*n])
        for j in eachindex(zq)
            total[c,j]=_barycentric_interpolate(response.z,values,bw,zq[j])
        end
    end
    incident=zeros(ComplexF64,4,length(zq))
    for j in eachindex(zq)
        phase=exp(im*Float64(kz)*(zq[j]-zedge))
        incident[:,j].=edge.*phase
    end
    induced=total-incident
    norms=[sqrt(sum(wq.*abs2.(induced[c,:]))) for c in 1:4]
    total_norms=[sqrt(sum(wq.*abs2.(total[c,:]))) for c in 1:4]
    incident_norms=[sqrt(sum(wq.*abs2.(incident[c,:]))) for c in 1:4]
    input_norm=norm(edge[1:3]); input_norm>0 || error("外缘速度输入范数不能为零")
    return (G_vel=sqrt(sum(abs2,norms[1:3]))/input_norm,
            G_ref=sqrt(sum(abs2,norms[1:3]))/input_norm,
            G_u=norms[1]/input_norm,
            G_v=norms[2]/input_norm,
            G_w=norms[3]/input_norm,
            G_pressure=norms[4]/input_norm,
            G_total_vel=sqrt(sum(abs2,total_norms[1:3]))/input_norm,
            G_D=sqrt(sum(abs2,total_norms[1:3]))/input_norm,
            G_incident_vel=sqrt(sum(abs2,incident_norms[1:3]))/input_norm,
            induced_norms=norms,
            total_norms=total_norms,
            incident_norms=incident_norms,
            input_norm=input_norm,
            z=zq,weights=wq,total=total,incident=incident,induced=induced,
            zi=zi,Nq=Nq,kz=Float64(kz))
end

"""用 companion 广义特征值问题求解二次空间模式。"""
function solve_spatial(op; shift=op.params.spatial_shift, nev=op.params.nev)
    n=size(op.L0,1); Id=Matrix{ComplexF64}(I,n,n); Z=zeros(ComplexF64,n,n)
    # y=[alpha*q; q]，因此第一行使用 [-L1 -L0]。
    A=[-op.L1 -op.L0; Id Z]; B=[op.L2 Z; Z Id]
    E=eigen(A,B); vals=ComplexF64.(E.values); vecs=ComplexF64.(E.vectors)
    ids=findall(i->isfinite(real(vals[i]))&&isfinite(imag(vals[i]))&&abs(vals[i])<10,eachindex(vals))
    sort!(ids,by=i->abs(vals[i]-shift))
    ids=ids[1:min(nev,length(ids))]
    return [(alpha=vals[i],vector=vecs[n+1:2n,i]) for i in ids]
end

"""
组装空间二次铅笔的伴随算子。

`formulation=:discrete` 是默认选项：对已经施加同一组边界约束的离散
二次铅笔逐项取转置。这样得到的是 notebook 使用的双线性离散伴随，
并且在相同离散系统上严格保持 LST/ALST 特征值一致。

`formulation=:formal` 按 notebook 的连续积分分部公式，从 `A,B,C,dC`
等块重新组装形式伴随。它适合核对连续推导；在半无限映射、非对称
Chebyshev 权重或不同边界删行下，形式伴随的离散谱不必与原离散谱
逐位相同。
"""
function assemble_adjoint_spatial(op; formulation::Symbol=:discrete)
    formulation in (:discrete,:formal) ||
        error("formulation 必须是 :discrete 或 :formal")
    if formulation == :discrete
        return (L0=transpose(op.L0), L1=transpose(op.L1), L2=transpose(op.L2),
                z=op.z, D=op.D, params=op.params, blocks=op.blocks,
                formulation=formulation)
    end

    b=op.blocks; p=op.params; n=p.N+1
    Kd=kron(Matrix{ComplexF64}(I,4,4),b.D)
    Kd2=kron(Matrix{ComplexF64}(I,4,4),b.D2)
    # 这是 receptivity.ipynb 中 assemble_adjmat 的四变量版本。
    A0=transpose(b.D1)+im*p.beta*p.R*transpose(b.B)-
       im*p.omega*transpose(b.Ta)-p.beta^2*p.R^2*transpose(b.Vyy)-
       transpose(b.dC)-im*p.beta*p.R*transpose(b.dVyz)+transpose(b.d2Vzz)-
       (transpose(b.C)+im*p.beta*p.R*transpose(b.Vyz)-
        2*transpose(b.dVzz))*Kd+transpose(b.Vzz)*Kd2
    A1=im*transpose(b.A)-p.beta*p.R*transpose(b.Vxy)-
       im*transpose(b.dVxz)-im*transpose(b.Vxz)*Kd
    A2=-transpose(b.Vxx)
    keep=b.keep
    return (L0=Matrix{ComplexF64}(A0[keep,keep]),
            L1=Matrix{ComplexF64}(A1[keep,keep]),
            L2=Matrix{ComplexF64}(A2[keep,keep]),
            z=b.z,D=b.D,params=p,blocks=b,formulation=formulation)
end

"""直接从基本流组装伴随空间算子。"""
function assemble_adjoint_spatial(base::BaseFlow;
                                  p::LSTParameters=LSTParameters(mode=:spatial),
                                  formulation::Symbol=:discrete)
    return assemble_adjoint_spatial(assemble_spatial(base;p=p);
                                    formulation=formulation)
end

"""求解伴随空间二次特征值问题。默认使用离散双线性转置。"""
function solve_adjoint_spatial(op; shift=op.params.spatial_shift,
                               nev=op.params.nev,
                               formulation::Symbol=:discrete)
    aop=assemble_adjoint_spatial(op; formulation=formulation)
    modes=solve_spatial(aop; shift=shift, nev=nev)
    return modes
end

"""从给定基本流直接求解伴随空间模式。"""
function solve_adjoint_spatial(base::BaseFlow;
                               p::LSTParameters=LSTParameters(mode=:spatial),
                               shift=p.spatial_shift, nev=p.nev,
                               formulation::Symbol=:discrete)
    op=assemble_spatial(base;p=p)
    return solve_adjoint_spatial(op; shift=shift, nev=nev,
                                 formulation=formulation)
end

"""从内部 von Kármán ODE 基流直接求解伴随空间模式。"""
function solve_adjoint_spatial(; p::LSTParameters=LSTParameters(mode=:spatial),
                               base_params::VonKarmanBaseParameters=VonKarmanBaseParameters(),
                               shift=p.spatial_shift, nev=p.nev,
                               formulation::Symbol=:discrete)
    base=solve_von_karman_baseflow(params=base_params)
    return solve_adjoint_spatial(base; p=p, shift=shift, nev=nev,
                                 formulation=formulation)
end

"""ALST 的简写接口，默认返回与离散 LST 同谱的双线性伴随。"""
assemble_alst(args...; kwargs...) = assemble_adjoint_spatial(args...; kwargs...)
solve_alst(args...; kwargs...) = solve_adjoint_spatial(args...; kwargs...)

"""计算原始变量二次铅笔的伴随模态归一化因子。

若 `weight` 给出，则使用 notebook 中的双线性权重
`transpose(q_adj) * weight * P'(alpha) * q`；否则使用单位权重。
"""
function adjoint_overlap(op, mode, adjoint_mode; weight=nothing)
    α=mode.alpha; αa=adjoint_mode.alpha
    derivative=op.L1+(α+αa)*op.L2
    left=transpose(adjoint_mode.vector)
    weight === nothing && return (left*derivative*mode.vector)[1]
    return (left*weight*derivative*mode.vector)[1]
end

"""用同一空间铅笔比较 LST 与离散 ALST 的目标特征值。"""
function compare_lsts_alst(op; shift=op.params.spatial_shift, nev=1)
    direct=solve_spatial(op; shift=shift, nev=nev)
    adjoint=solve_adjoint_spatial(op; shift=shift, nev=nev,
                                  formulation=:discrete)
    isempty(direct) && error("LST 没有找到有限空间特征值")
    isempty(adjoint) && error("ALST 没有找到有限空间特征值")
    return (alpha=direct[1].alpha,
            alpha_adjoint=adjoint[1].alpha,
            difference=adjoint[1].alpha-direct[1].alpha,
            overlap=adjoint_overlap(op,direct[1],adjoint[1]),
            lst=direct[1], alst=adjoint[1])
end

function solve_lsts(base::BaseFlow; p::LSTParameters=LSTParameters())
    p.mode==:temporal && return solve_temporal(assemble_temporal(base;p=p))
    p.mode==:spatial && return solve_spatial(assemble_spatial(base;p=p))
    error("mode 必须是 :temporal 或 :spatial")
end

"""直接从 von Kármán ODE 基流开始求解，不需要外部基流文件。"""
function solve_lsts(; p::LSTParameters=LSTParameters(),
                    base_params::VonKarmanBaseParameters=VonKarmanBaseParameters())
    base = solve_von_karman_baseflow(; params=base_params)
    return solve_lsts(base; p=p)
end

end # module VonKarmanLST

# ============================================================================
# 使用说明
# ============================================================================
#
# 1. 加载求解器
#
#    include("VonKarmanLST.jl")
#    using .VonKarmanLST
#
# 2. 求解 von Kármán 相似性基本流
#
#    base_parameters = VonKarmanBaseParameters(
#        eta_max=30.0,          # ODE 积分终点；空间半无限映射建议至少取 20
#        n_points=6001,         # 基流保存点数
#        initial_fprime=0.51,   # F'(0) 的 Newton 初值
#        initial_gprime=-0.616, # G'(0) 的 Newton 初值
#        newton_tolerance=1e-11,
#        newton_max_iterations=20,
#    )
#    base = solve_von_karman_baseflow(params=base_parameters)
#
#    返回的 BaseFlow 使用随圆盘旋转坐标系：
#        U = F,
#        V = G - 1,
#        W = H。
#
# 3. 时间模式
#
#    固定 alpha、beta 和 R，求时间频率 omega：
#
#    p_t = LSTParameters(
#        mode=:temporal,
#        R=285.36,
#        alpha=0.38482,
#        beta=0.07759,
#        Ro=1.0,
#        Co=2.0,
#        N=99,
#        temporal_zmax=50.0,
#    )
#
#    op_t = assemble_temporal(base; p=p_t)
#    modes_t = solve_temporal(op_t)
#    omega_1 = modes_t[1].omega
#
# 4. 空间模式
#
#    固定 omega、beta 和 R，求径向波数 alpha：
#
#    p_s = LSTParameters(
#        mode=:spatial,
#        R=285.36,
#        omega=0.0,
#        beta=0.07759,
#        Ro=1.0,
#        Co=2.0,
#        N=99,
#        spatial_shift=0.38 + 0im,
#        nev=1,
#    )
#
#    op_s = assemble_spatial(base; p=p_s)
#    modes_s = solve_spatial(op_s)
#    alpha_1 = modes_s[1].alpha
#
#    空间模式使用 z∈[0,+∞) 的半无限 Chebyshev 映射，不需要 zmax。
#
# 5. 伴随空间模式
#
#    op_s = assemble_spatial(base; p=p_s)
#    modes_a = solve_alst(op_s)
#    check = compare_lsts_alst(op_s)
#
#    默认 ALST 是施加同一边界约束后的离散双线性转置，因此
#    check.alpha 和 check.alpha_adjoint 应在机器精度内一致。
#    notebook 风格的连续形式伴随可用：
#
#    op_a_formal = assemble_alst(op_s; formulation=:formal)
#
# 6. 光滑盘面直接受迫响应
#
#    外缘输入必须是满足外部连续性/传播关系的复幅值：
#        edge_state = ComplexF64[u_e, v_e, w_e, p_e]
#    这里的 alpha=p.alpha 是外部给定波数，不是重新求的 LST 特征值。
#
#    p_f = LSTParameters(
#        mode=:spatial, R=285.36, alpha=0.38482,
#        beta=0.07759, omega=0.0, N=49,
#    )
#    response = solve_forced_spatial(
#        base, edge_state; p=p_f, zmax=20.0, lift=:smooth,
#    )
#    total_q = response.q
#    norm(response.interior_residual)  # 应接近离散误差
#    norm(response.boundary_error)     # 应接近离散误差
#
#    `lift=:linear` 和 `lift=:smooth` 应给出相同的重构总场；
#    只有提升场和齐次修正场分别依赖提升函数。
#
# 7. 自由涡波极化和响应增益
#
#    k = [alpha_e, beta_e, k_z]
#    states = vortical_edge_states(k; amplitude=1e-3)
#    response = solve_forced_spatial(base, states[1]; p=p_f, zmax=30.0)
#    gains = response_gains(response; edge_state=states[1], zi=10.0, Nq=128)
#
#    `gains.G_vel` 是扣除物理入射场后的诱导速度增益，
#    `gains.G_total_vel`/`gains.G_D` 是总场速度增益；两者必须分开报告。
#    `gains.G_u,G_v,G_w` 是诱导速度分量增益。
#
#    远场一致性检查：
#        modes_inf = solve_farfield_modes(base; p=p_f)
#        M0 = farfield_matrix(base; p=p_f, kz=0.0+0im)
#        diagnostics = farfield_mode_diagnostics(base; p=p_f)
#    `diagnostics` 按波动涡量比例标记 potential_like/vortical_like。
#
# 8. 统一接口
#
#    时间模式：
#        modes_t = solve_lsts(base; p=p_t)
#
#    空间模式：
#        modes_s = solve_lsts(base; p=p_s)
#
#    也可以让统一接口内部自动求解基本流：
#
#        modes_t = solve_lsts(p=p_t, base_params=base_parameters)
#        modes_s = solve_lsts(p=p_s, base_params=base_parameters)
#
# 9. 文献中性点参数
#
#    Malik Type-I 上支中性点：
#        R       = 285.36
#        alpha   = 0.38482
#        beta    = 0.07759
#        omega   = 0
#        Ro      = 1
#        Co      = 2
#
#    以默认 eta_max=30、原始变量半无限配置计算时，Type-I 空间结果应约为
#        alpha = 0.38505249 + 0.00034867im
#    这是当前离散原始变量算子的收敛值，不应宣称严格等于 Malik 的
#    0.38482+0im；空间残差和网格收敛只能说明本算子求解稳定。
#
# 10. 注意事项
#
#    * 时间模式的特征值是 omega。
#    * 空间模式的特征值是 alpha。
#    * 两种模式共用同一个内部 ODE 基流，但离散矩阵算子不同。
#    * 本文件不读取外部基流文件，也不自动写入 CSV、DAT 或报告文件。

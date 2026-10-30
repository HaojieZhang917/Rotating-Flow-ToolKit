#=
六维一阶 frozen-R LNS 求解器。

状态选为 Y=[u,v,w,p,u_z,v_z]。本文件从 VonKarmanLST.jl 已有的四变量
原始算子块逐点构造 Y_z=B(z)Y，不改变旧四变量求解器。
=#
module SixStateLNS

using LinearAlgebra
using ..VonKarmanLST

export assemble_six_state_lns, solve_six_state_lns,
       farfield_six_state_modes, modal_outer_boundary,
       six_state_residuals, state_to_primitive,
       as_four_state_response, solve_six_state_forced,
       solve_six_state_coupled_response, six_state_forcing

"从全局 4n×4n 块矩阵提取第 j 个节点的 4×4 局部块。"
function _local_block(M,j,n)
    rows=[(r-1)*n+j for r in 1:4]
    return ComplexF64.(M[rows,rows])
end

"构造一个节点的六维一阶矩阵，严格使用四变量原始算子。"
function _local_six_state_matrix(blocks,p,j)
    n=p.N+1
    A=_local_block(blocks.A,j,n)
    B=_local_block(blocks.B,j,n)
    C=_local_block(blocks.C,j,n)
    Ta=_local_block(blocks.Ta,j,n)
    D1=_local_block(blocks.D1,j,n)
    Vxx=_local_block(blocks.Vxx,j,n)
    Vyy=_local_block(blocks.Vyy,j,n)
    Vzz=_local_block(blocks.Vzz,j,n)
    α=p.alpha; β=p.beta; R=p.R
    K=D1+im*β*R*B-im*p.omega*Ta+im*α*A-
      α^2*Vxx-β^2*R^2*Vyy

    # q_z 的六维表示。连续性只含 w_z，因此可直接消去 w_z。
    abs(C[4,3])>1e-13 || error("连续性方程的 w_z 系数退化")
    g=zeros(ComplexF64,6)
    g[1:4].=-K[4,:]./C[4,3]
    abs(g[3])<1e-12 && abs(g[4])<1e-12 ||
        error("当前连续性消元出现 w 或 p 的隐式导数，六维显式形式不适用")

    qmap=zeros(ComplexF64,4,6)
    qmap[1:4,1:4].=Matrix{ComplexF64}(I,4,4)
    qzmap=zeros(ComplexF64,4,6)
    qzmap[1,5]=1
    qzmap[2,6]=1
    qzmap[3,:]=g
    # w_zz 由连续性对 z 求导得到。
    wzz=zeros(ComplexF64,6)
    wzz[5]=g[1]
    wzz[6]=g[2]

    # 第三个动量方程恢复 p_z。
    abs(C[3,4])>1e-13 || error("法向动量方程的 p_z 系数退化")
    qzmap[4,:]=vec(-(Vzz[3,3].*reshape(wzz,1,6) .+
                     reshape(C[3,1:3],1,3)*qzmap[1:3,:] .+
                     reshape(K[3,:],1,4)*qmap)./C[3,4])

    # 前两个动量方程恢复 u_zz、v_zz。
    abs(Vzz[1,1])>1e-13 && abs(Vzz[2,2])>1e-13 ||
        error("切向动量方程的二阶法向导数系数退化")
    uzz=vec(-(reshape(C[1,:],1,4)*qzmap+reshape(K[1,:],1,4)*qmap)./Vzz[1,1])
    vzz=vec(-(reshape(C[2,:],1,4)*qzmap+reshape(K[2,:],1,4)*qmap)./Vzz[2,2])

    B6=zeros(ComplexF64,6,6)
    B6[1,:].=qzmap[1,:]
    B6[2,:].=qzmap[2,:]
    B6[3,:].=qzmap[3,:]
    B6[4,:].=qzmap[4,:]
    B6[5,:].=uzz
    B6[6,:].=vzz
    return B6,K,C,Vzz
end

"""
组装六维 Chebyshev 配点系统。

返回的 `A_dynamic` 对应 `Y_z-B(z)Y=0`，节点变量按六个状态分块排列。
"""
function assemble_six_state_lns(base::VonKarmanLST.BaseFlow;
                                p::VonKarmanLST.LSTParameters=VonKarmanLST.LSTParameters(mode=:spatial),
                                zmax::Float64=30.0)
    blocks=VonKarmanLST.assemble_spatial_blocks(base;p=p,grid=:affine,zmax=zmax)
    n=p.N+1
    D=blocks.D
    Bnodes=Vector{Matrix{ComplexF64}}(undef,n)
    Kblocks=Vector{Matrix{ComplexF64}}(undef,n)
    Cblocks=Vector{Matrix{ComplexF64}}(undef,n)
    Vblocks=Vector{Matrix{ComplexF64}}(undef,n)
    G=zeros(ComplexF64,6n,6n)
    for j in 1:n
        B6,K,C,Vzz=_local_six_state_matrix(blocks,p,j)
        Bnodes[j]=B6; Kblocks[j]=K; Cblocks[j]=C; Vblocks[j]=Vzz
        for r in 1:6, c in 1:6
            G[(r-1)*n+j,(c-1)*n+j]=B6[r,c]
        end
    end
    D6=kron(Matrix{ComplexF64}(I,6,6),D)
    A_dynamic=D6-G
    return (A_dynamic=A_dynamic,D=D,z=blocks.z,Bnodes=Bnodes,
            Kblocks=Kblocks,Cblocks=Cblocks,Vzzblocks=Vblocks,
            blocks=blocks,params=p,zmax=zmax)
end

"将六维状态向量恢复为节点上的四变量原始场。"
function state_to_primitive(Y::AbstractVector,n::Int)
    length(Y)==6n || error("六维状态向量长度必须为 6n")
    return ComplexF64.(Y[1:4n])
end

"""
把六维响应转换为旧二阶耦合接口使用的四变量响应对象。

`qz` 中的法向导数优先来自六维状态方程；同时返回 `qz_cheb` 供验证
与 Chebyshev 求导结果的差异。这样不会重复引入外部幅值或粗糙度幅值。
"""
function as_four_state_response(sol)
    n=length(sol.z); Y=ComplexF64.(sol.Y)
    q=Y[1:4n]
    qz=zeros(ComplexF64,4n)
    qz[1:n].=Y[4n+1:5n]
    qz[n+1:2n].=Y[5n+1:6n]
    for j in 1:n
        yd=sol.operator.Bnodes[j]*ComplexF64[Y[(c-1)*n+j] for c in 1:6]
        qz[2n+j]=yd[3]
        qz[3n+j]=yd[4]
    end
    qz_cheb=vcat(sol.D*q[1:n],sol.D*q[n+1:2n],sol.D*q[2n+1:3n],sol.D*q[3n+1:4n])
    shear=ComplexF64[qz[1],qz[n+1],qz[2n+1]]
    return (q=q,qz=qz,qz_explicit=qz,z=sol.z,D=sol.D,params=sol.params,
            zmax=sol.zmax,wall_shear=shear,qz_cheb=qz_cheb,
            velocity_derivative_difference=norm(qz[1:3n]-qz_cheb[1:3n],Inf),
            derivative_difference=norm(qz-qz_cheb,Inf),
            source=zeros(ComplexF64,4n),state_solution=sol)
end

"计算六维方程残差和恢复后的四变量 LNS 残差。"
function six_state_residuals(op,Y::AbstractVector; forcing=zeros(ComplexF64,length(Y)))
    n=length(op.z)
    length(Y)==6n || error("状态向量长度与算子不一致")
    length(forcing)==length(Y) || error("forcing 长度与状态不一致")
    dynamic=op.A_dynamic*Y-ComplexF64.(forcing)
    primitive=state_to_primitive(Y,n)
    # 用六维导数恢复 q_z 和 q_zz，逐点计算原始四变量方程。
    qres=zeros(ComplexF64,4n)
    for j in 1:n
        q=primitive[[(r-1)*n+j for r in 1:4]]
        # 直接从 Y_z=B6Y 读取 uz,vz,wz,pz。
        yd=op.Bnodes[j]*(ComplexF64[Y[(r-1)*n+j] for r in 1:6])
        qz=ComplexF64[ Y[4n+j],Y[5n+j],yd[3],yd[4] ]
        qzz=ComplexF64[yd[5],yd[6],
                      op.Bnodes[j][3,1]*Y[4n+j]+op.Bnodes[j][3,2]*Y[5n+j],0]
        rloc=op.Vzzblocks[j]*qzz+op.Cblocks[j]*qz+op.Kblocks[j]*q
        for r in 1:4; qres[(r-1)*n+j]=rloc[r]; end
    end
    return (dynamic=dynamic,primitive=qres,
            dynamic_inf=norm(dynamic,Inf),primitive_inf=norm(qres,Inf))
end

"从远场六个有限模态中提取与六维状态一致的 [q;u_z;v_z] 状态。"
function farfield_six_state_modes(base::VonKarmanLST.BaseFlow;
                                  p::VonKarmanLST.LSTParameters=VonKarmanLST.LSTParameters(mode=:spatial),
                                  max_abs_kz::Real=100.0)
    full=VonKarmanLST.solve_farfield_state_modes(base;p=p,max_abs_kz=max_abs_kz)
    modes=NamedTuple[]
    for m in full.modes
        y6=ComplexF64.(vcat(m.q,m.qz[1:2]))
        push!(modes,(kz=m.kz,lambda=m.lambda,state=y6,
                     q=m.q,qz=m.qz,pencil_residual=m.pencil_residual,
                     quadratic_residual=m.quadratic_residual))
    end
    return modes
end

"用 SVD 构造三个有界有限模态的左零空间边界矩阵。"
function modal_outer_boundary(base::VonKarmanLST.BaseFlow;
                              p::VonKarmanLST.LSTParameters=VonKarmanLST.LSTParameters(mode=:spatial),
                              imag_tolerance::Real=1e-9)
    modes=farfield_six_state_modes(base;p=p)
    decaying=[m for m in modes if imag(m.kz)>imag_tolerance]
    length(decaying)==3 || error("当前参数下有界有限模态数为 $(length(decaying))，不是三维")
    V=hcat([m.state for m in decaying]...)
    F=svd(V;full=true)
    C=adjoint(F.U[:,4:6])
    growing=[m for m in modes if imag(m.kz)<-imag_tolerance]
    Vg=hcat([m.state for m in growing]...)
    return (C=C,V=V,modes=decaying,growing_modes=growing,
            rank=rank(C),growth_constraint_rank=rank(C*Vg),
            constraint_residual=norm(C*V),condition=cond(V),all_modes=modes)
end

function _incident_state(base,p; z::Float64, zref::Float64, amplitude::Complex,
                         kz_target::Complex=-0.1491006263im)
    modes=farfield_six_state_modes(base;p=p)
    idx=argmin([abs(m.kz-kz_target) for m in modes])
    m=modes[idx]
    y=copy(m.state)
    vnorm=norm(y[1:3]); vnorm>0 || error("入射模态速度幅值为零")
    y./=vnorm
    return amplitude.*y.*exp(im*m.kz*(z-zref)),m
end

"""
求解六维一阶受迫问题。

`boundary=:dirichlet` 在外缘施加三分量速度 Dirichlet 条件；
`boundary=:modal` 施加衰减子空间的三条模态条件。两者都只使用
同样的壁面三速度条件，不指定壁面或外缘压力。
"""
function solve_six_state_lns(base::VonKarmanLST.BaseFlow;
                             p::VonKarmanLST.LSTParameters=VonKarmanLST.LSTParameters(mode=:spatial),
                             zmax::Float64=30.0,
                             boundary::Symbol=:dirichlet,
                             incident_state=nothing,
                             zref::Float64=10.0,
                             incident_amplitude::Complex=1e-3+0im,
                             edge_velocity=nothing)
    boundary in (:dirichlet,:modal) || error("boundary 必须是 :dirichlet 或 :modal")
    op=assemble_six_state_lns(base;p=p,zmax=zmax)
    n=p.N+1; L=op.z[end]
    yinc,m_inc=incident_state===nothing ?
        _incident_state(base,p;z=L,zref=zref,amplitude=incident_amplitude) :
        (ComplexF64.(incident_state),nothing)
    length(yinc)==6 || error("incident_state 必须为六维外缘状态")
    A=copy(op.A_dynamic); rhs=zeros(ComplexF64,6n)
    wall_rows=[1,n+1,2n+1]
    outer_rows=[n,2n,3n]
    # 壁面 u=v=w=0。
    for (r,c) in zip(wall_rows,[1,n+1,2n+1])
        A[r,:].=0; A[r,c]=1; rhs[r]=0
    end
    bc_data=nothing
    if boundary==:dirichlet
        edge=edge_velocity===nothing ? yinc[1:3] : ComplexF64.(edge_velocity)
        length(edge)==3 || error("edge_velocity 必须包含 [u,v,w]")
        for (r,c,val) in zip(outer_rows,[n,2n,3n],edge)
            A[r,:].=0; A[r,c]=1; rhs[r]=val
        end
        bc_data=(kind=:dirichlet,edge_velocity=edge)
    else
        mb=modal_outer_boundary(base;p=p)
        C=mb.C
        for k in 1:3
            r=outer_rows[k]; A[r,:].=0
            for c in 1:6; A[r,(c-1)*n+n]=C[k,c]; end
            # C 已经是左零空间矩阵；边界方程使用普通矩阵乘法，
            # 不能再用 dot（dot 会额外共轭 C）。
            rhs[r]=sum(C[k,:].*yinc)
        end
        bc_data=(kind=:modal,C=C,modal=mb,scattering_state=zeros(ComplexF64,6))
    end
    Y=A\rhs
    outer_state=ComplexF64[Y[n],Y[2n],Y[3n],Y[4n],Y[5n],Y[6n]]
    if boundary==:modal
        bc_data=(kind=:modal,C=bc_data.C,modal=bc_data.modal,
                 scattering_state=outer_state-yinc)
    end
    return (Y=Y,operator=op,A=A,rhs=rhs,z=op.z,D=op.D,
            params=p,zmax=zmax,boundary=boundary,
            incident_state=yinc,incident_mode=m_inc,
            outer_state=outer_state,bc_data=bc_data,
            residuals=six_state_residuals(op,Y))
end

"""
把原四变量方程右端 `[s_u,s_v,s_w,s_c]` 映射到六维一阶系统右端。

当前二阶非线性源项满足 `s_c=0`。连续性源若非零还需要其法向导数
才能构造 `w_zz`，因此本接口对非零连续性源明确报错，避免无声地遗漏项。
"""
function six_state_forcing(op,source::AbstractVector)
    n=length(op.z); length(source)==4n || error("source 必须包含 4n 个分量")
    s=ComplexF64.(source); sc=s[3n+1:4n]
    maximum(abs,sc)<=1e-12 || error("当前六维受迫接口要求连续性源为零；非零 s_c 需要额外的法向导数项")
    f=zeros(ComplexF64,6n)
    for j in 1:n
        C=op.Cblocks[j]; V=op.Vzzblocks[j]
        abs(C[3,4])>1e-13 || error("法向动量 p_z 系数退化")
        abs(V[1,1])>1e-13 && abs(V[2,2])>1e-13 || error("切向动量二阶导数系数退化")
        # 状态分块为 [u,v,w,p,u_z,v_z]；当前 s_c=0，前三个
        # 强迫分别进入 p_z、u_zz、v_zz。
        f[3n+j]=s[2n+j]/C[3,4]
        f[4n+j]=s[j]/V[1,1]
        f[5n+j]=s[n+j]/V[2,2]
    end
    return f
end

"""求解六维一阶受迫 LNS，支持三速度 Dirichlet 或 Modal 外缘条件。"""
function solve_six_state_forced(base::VonKarmanLST.BaseFlow,source::AbstractVector;
                                p::VonKarmanLST.LSTParameters=VonKarmanLST.LSTParameters(mode=:spatial),
                                zmax::Float64=30.0,
                                wall_velocity=zeros(ComplexF64,3),
                                boundary::Symbol=:modal,
                                edge_velocity=zeros(ComplexF64,3))
    boundary in (:modal,:dirichlet) || error("boundary 必须是 :modal 或 :dirichlet")
    op=assemble_six_state_lns(base;p=p,zmax=zmax); n=length(op.z)
    f=six_state_forcing(op,source); A=copy(op.A_dynamic); rhs=f
    wall=ComplexF64.(collect(wall_velocity)); length(wall)==3 || error("wall_velocity 必须有三个分量")
    edge=ComplexF64.(collect(edge_velocity)); length(edge)==3 || error("edge_velocity 必须有三个分量")
    wall_rows=[1,n+1,2n+1]; outer_rows=[n,2n,3n]
    for (r,c,val) in zip(wall_rows,[1,n+1,2n+1],wall)
        A[r,:].=0; A[r,c]=1; rhs[r]=val
    end
    modal=nothing
    if boundary==:dirichlet
        for (r,c,val) in zip(outer_rows,[n,2n,3n],edge)
            A[r,:].=0; A[r,c]=1; rhs[r]=val
        end
    else
        modal=modal_outer_boundary(base;p=p)
        for k in 1:3
            r=outer_rows[k]; A[r,:].=0
            for c in 1:6; A[r,(c-1)*n+n]=modal.C[k,c]; end
            rhs[r]=0
        end
    end
    Y=A\rhs
    outer=ComplexF64[Y[n],Y[2n],Y[3n],Y[4n],Y[5n],Y[6n]]
    return (Y=Y,operator=op,A=A,rhs=rhs,z=op.z,D=op.D,params=p,zmax=zmax,
            boundary=boundary,source=ComplexF64.(source),wall_velocity=wall,
            outer_state=outer,modal=modal,residuals=six_state_residuals(op,Y;forcing=f))
end

"""用六维受迫系统求解体积、壁面和完整二阶响应。"""
function solve_six_state_coupled_response(base::VonKarmanLST.BaseFlow,qh,qf;
                                          p_h=nothing,p_f=nothing,
                                          roughness_amplitude::Number=1.0+0im,
                                          zmax::Float64=30.0,
                                          boundary::Symbol=:modal)
    forcing=VonKarmanLST.nonlinear_cross_forcing(qh,qf,base;p_h=p_h,p_f=p_f)
    ph=forcing.params_h; pf=forcing.params_f
    pt=VonKarmanLST.LSTParameters(mode=:spatial,R=pf.R,alpha=forcing.target_alpha,
        beta=forcing.target_beta,omega=forcing.target_omega,Ro=pf.Ro,Co=pf.Co,
        N=pf.N,temporal_zmax=pf.temporal_zmax,spatial_shift=pf.spatial_shift,nev=pf.nev)
    scale=ComplexF64(roughness_amplitude); source=scale.*forcing.source
    shear=VonKarmanLST.wall_shear_response(qf); wall_velocity=-scale.*ComplexF64.(shear)
    zero_source=zeros(ComplexF64,length(source)); zero_wall=zeros(ComplexF64,3); zero_edge=zeros(ComplexF64,3)
    volume=solve_six_state_forced(base,source;p=pt,zmax=zmax,wall_velocity=zero_wall,boundary=boundary,edge_velocity=zero_edge)
    wall=solve_six_state_forced(base,zero_source;p=pt,zmax=zmax,wall_velocity=wall_velocity,boundary=boundary,edge_velocity=zero_edge)
    total=solve_six_state_forced(base,source;p=pt,zmax=zmax,wall_velocity=wall_velocity,boundary=boundary,edge_velocity=zero_edge)
    return (response=total,volume_response=volume,wall_response=wall,source=source,
            unit_source=forcing.source,wall_velocity=wall_velocity,unit_wall_velocity=-shear,
            forcing=forcing,params=pt,roughness_amplitude=scale)
end

end # module SixStateLNS

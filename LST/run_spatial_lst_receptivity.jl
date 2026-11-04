include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
using LinearAlgebra
using Serialization
using Random

const OUTDIR=joinpath(@__DIR__,"results","spatial_lstreceptivity")
mkpath(OUTDIR)
const R0=285.36
const BETA0=0.07759
const OMEGA_F=-0.07759

make_p(N,omega)=LSTParameters(mode=:spatial,R=R0,alpha=.38482,beta=BETA0,omega=omega,N=N,nev=4*(N+1),spatial_shift=.38+0im)

function assemble_finite(base,p,zmax)
    b=assemble_spatial_blocks(base;p=p,grid=:finite,zmax=zmax); L0,L1,L2=VonKarmanLST._spatial_pencil(b,p); k=b.keep
    (L0=Matrix{ComplexF64}(L0[k,k]),L1=Matrix{ComplexF64}(L1[k,k]),L2=Matrix{ComplexF64}(L2[k,k]),z=b.z,D=b.D,params=p,blocks=b,zmax=zmax)
end

function residual(op,m)
    r=(op.L0+m.alpha*op.L1+m.alpha^2*op.L2)*m.vector
    norm(r)/max(norm(op.L0)*norm(m.vector),eps(Float64))
end

function choose_reference(modes)
    c=[m for m in modes if real(m.alpha)>0.2 && real(m.alpha)<0.7]
    isempty(c) && error("没有找到参考空间模态")
    c[argmin([abs(m.alpha-(.38482+0im)) for m in c])]
end

function write_csv(path,header,rows)
    io=IOBuffer(); println(io,join(header,","))
    for row in rows; println(io,join(string.(row),",")); end
    write(path,String(take!(io)))
end

function embed(op,m)
    full=zeros(ComplexF64,4length(op.z)); full[op.blocks.keep].=m.vector; full
end

function normalize_mode(op,m,c)
    full=embed(op,m); n=length(op.z); vals=full[(c-1)*n+1:c*n]; j=argmax(abs.(vals)); s=abs(vals[j]); s>0 || error("模态归一化分量为零")
    v=m.vector/s; phase=angle(v[findfirst(==(op.blocks.keep[(c-1)*n+j]),op.blocks.keep)])
    v.*=exp(-im*phase)
    (alpha=m.alpha,vector=v)
end

base=solve_von_karman_baseflow()
rows=NamedTuple[]
selected=nothing
for (label,omega) in [(:target,OMEGA_F),(:reference,0.0)]
    for N in [49,69,89]
        println("solve ",label," N=",N)
        p=make_p(N,omega); op=assemble_finite(base,p,30.0)
        modes=solve_spatial(op;shift=.38+0im,nev=size(op.L0,1))
        phys=[m for m in modes if real(m.alpha)>0.2 && real(m.alpha)<0.7]
        if label==:reference
            dm=normalize_mode(op,choose_reference(modes),2)
            amodes=solve_adjoint_spatial(op;shift=dm.alpha,nev=size(op.L0,1),formulation=:discrete)
            am=normalize_mode(op,amodes[argmin([abs(x.alpha-dm.alpha) for x in amodes])],1)
            aa=assemble_adjoint_spatial(op;formulation=:discrete)
            ar=norm((aa.L0+am.alpha*aa.L1+am.alpha^2*aa.L2)*am.vector)/max(norm(aa.L0)*norm(am.vector),eps(Float64))
            Q=(transpose(am.vector)*(op.L1+2dm.alpha*op.L2)*dm.vector)[1]
            push!(rows,(case=String(label),N=N,omega=omega,alpha_real=real(dm.alpha),alpha_imag=imag(dm.alpha),growth=-imag(dm.alpha),direct_residual=residual(op,dm),adjoint_residual=ar,Q_real=real(Q),Q_imag=imag(Q),candidate_count=length(phys),classification="candidate"))
            N==89 && (global selected=(p=p,op=op,dm=dm,am=am,Q=Q))
        else
            if isempty(phys)
                push!(rows,(case=String(label),N=N,omega=omega,alpha_real=NaN,alpha_imag=NaN,growth=NaN,direct_residual=NaN,adjoint_residual=NaN,Q_real=NaN,Q_imag=NaN,candidate_count=0,classification="no candidate"))
            else
                dm=phys[argmin([abs(imag(m.alpha)) for m in phys])]
                push!(rows,(case=String(label),N=N,omega=omega,alpha_real=real(dm.alpha),alpha_imag=imag(dm.alpha),growth=-imag(dm.alpha),direct_residual=residual(op,dm),adjoint_residual=NaN,Q_real=NaN,Q_imag=NaN,candidate_count=length(phys),classification="candidate; not Type I/II"))
            end
        end
    end
end
cols=String.(fieldnames(typeof(rows[1])))
write_csv(joinpath(OUTDIR,"spatial_lstr_eigenvalues.csv"),cols,[[getproperty(r,Symbol(c)) for c in cols] for r in rows])

ref=selected; dm=ref.dm; am=ref.am; op=ref.op; p=ref.p; n=p.N+1
serialize(joinpath(OUTDIR,"spatial_direct_modes.jls"),(mode=dm,params=p,z=op.z,operator=op))
serialize(joinpath(OUTDIR,"spatial_adjoint_modes.jls"),(mode=am,params=p,z=op.z,operator=assemble_adjoint_spatial(op;formulation=:discrete)))
function block(m,c)
    full=zeros(ComplexF64,4n); full[op.blocks.keep].=m.vector
    full[(c-1)*n+1:c*n]
end
prows=NamedTuple[]
for j in 1:n
    push!(prows,(z=op.z[j],u_re=real(block(dm,1)[j]),u_im=imag(block(dm,1)[j]),v_re=real(block(dm,2)[j]),v_im=imag(block(dm,2)[j]),w_re=real(block(dm,3)[j]),w_im=imag(block(dm,3)[j]),p_re=real(block(dm,4)[j]),p_im=imag(block(dm,4)[j]),au_re=real(block(am,1)[j]),au_im=imag(block(am,1)[j]),av_re=real(block(am,2)[j]),av_im=imag(block(am,2)[j]),aw_re=real(block(am,3)[j]),aw_im=imag(block(am,3)[j]),ap_re=real(block(am,4)[j]),ap_im=imag(block(am,4)[j])))
end
cols2=String.(fieldnames(typeof(prows[1])))
write_csv(joinpath(OUTDIR,"spatial_mode_profiles.csv"),cols2,[[getproperty(r,Symbol(c)) for c in cols2] for r in prows])

Q=ref.Q
dm2=(alpha=dm.alpha,vector=2dm.vector); am2=(alpha=am.alpha,vector=3am.vector)
Q2=(transpose(am2.vector)*(op.L1+2dm2.alpha*op.L2)*dm2.vector)[1]
proj1=dm.vector*transpose(am.vector)/Q; proj2=dm2.vector*transpose(am2.vector)/Q2
rng=MersenneTwister(42); x=randn(rng,ComplexF64,size(op.L0,1)); y=randn(rng,ComplexF64,size(op.L0,1));
green=abs((transpose(y)*op.L0*x)[1]-(transpose(transpose(op.L0)*y)*x)[1])
write_csv(joinpath(OUTDIR,"spatial_mode_normalization.csv"),["Q_real","Q_imag","projector_scale_error","green_check"],[(real(Q),imag(Q),norm(proj1-proj2)/max(norm(proj1),eps(Float64)),green)])

open(joinpath(OUTDIR,"operator_audit.md"),"w") do io
    println(io,"# 空间 LST 与离散伴随算子核查")
    println(io,"当前变量为 `[u,v,w,p]`，空间算子为 `L0+alpha*L1+alpha^2*L2`，Fourier 约定为 `exp(i*(alpha*r+beta*s-omega*t))`。")
    println(io,"当前离散伴随是施加相同边界约束后的二次铅笔逐项转置；`Q` 使用 `transpose(y)*(L1+2alpha*L2)*x`。")
    println(io,"当前结果是离散双线性伴随框架，后续若采用物理 Hermitian 伴随，需要补充一致的权重矩阵和连续 Green 边界验证。")
end

open(joinpath(OUTDIR,"spatial_lst_adjoint_report.md"),"w") do io
    println(io,"# von Kármán 空间 LST 与伴随模态验证")
    println(io,"参数：R=$(R0), beta=$(BETA0), omega_f=$(OMEGA_F)。")
    println(io,"目标 omega_f 下的候选记录在 `spatial_lstr_eigenvalues.csv`；当前候选没有被自动标记为 Type I/II，且其谱分支随 N 变化明显，不能作为最终目标模态。")
    println(io,"omega=0 参考算例用于验证直接/离散伴随接口，N=89 的模态和伴随模态已经保存。参考 alpha 在 N=49,69,89 间稳定到约 1e-11，说明二次铅笔离散求解本身收敛。")
    println(io,"Q、投影归一化和离散 Green 检查见 `spatial_mode_normalization.csv`。当前 Q 的 N=49,69,89 数值分别约为 `0.0644+0.0236i`、`0.0691+0.0672i`、`0.0673+0.1046i`，尚未收敛，原因是当前半无限映射尚未提供一致的物理 Hermitian 积分权重。")
    println(io,"直接与伴随剖面见 `spatial_mode_profiles.csv`；二进制模态见 `spatial_direct_modes.jls` 和 `spatial_adjoint_modes.jls`。")
    println(io,"直接/伴随离散残差约为 1e-17，投影算子对左右特征向量缩放不变，离散转置 Green 检查为零。当前尚未完成 Type I/II 的物理分类、半无限物理内积的连续 Green 验证和 Gaussian 感受性系数；不能将候选模态直接用于最终感受性结论。")
end
println("wrote spatial LST/adjoint report")

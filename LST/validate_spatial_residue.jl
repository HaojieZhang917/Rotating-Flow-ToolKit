include(joinpath(@__DIR__,"VonKarmanLST.jl")); using .VonKarmanLST
using LinearAlgebra, Serialization

const OUTDIR=joinpath(@__DIR__,"results","spatial_adjoint_validation"); mkpath(OUTDIR)
const R0=285.36; const BETA=.07759; const OMEGA=0.0; const ZMAX=30.0
make_p(N)=LSTParameters(mode=:spatial,R=R0,alpha=.38482,beta=BETA,omega=OMEGA,N=N,nev=4*(N+1),spatial_shift=.38+0im)

function make_op(base,p)
    b=assemble_spatial_blocks(base;p=p,grid=:finite,zmax=ZMAX); a0,a1,a2=VonKarmanLST._spatial_pencil(b,p); k=b.keep
    (L0=Matrix{ComplexF64}(a0[k,k]),L1=Matrix{ComplexF64}(a1[k,k]),L2=Matrix{ComplexF64}(a2[k,k]),rawL0=a0,rawL1=a1,rawL2=a2,keep=k,z=b.z,D=b.D,params=p,blocks=b)
end
function embed(op,m)
    x=zeros(ComplexF64,4length(op.z)); x[op.keep].=m.vector; x
end
function normalize(op,m,c)
    n=length(op.z); f=embed(op,m); a=f[(c-1)*n+1:c*n]; j=argmax(abs.(a)); s=abs(a[j]); x=m.vector/s; full=zeros(ComplexF64,4n); full[op.keep].=x; phase=angle(full[(c-1)*n+j]); x.*=exp(-im*phase); (alpha=m.alpha,vector=x)
end
function choose(ms)
    c=[m for m in ms if real(m.alpha)>.2&&real(m.alpha)<.7]; c[argmin([abs(m.alpha-.38482) for m in c])]
end
function err_for(op,x,y,Q,b,δ)
    M=op.L0+(x.alpha+δ)*op.L1+(x.alpha+δ)^2*op.L2
    q=δ.*(M\b); pred=x.vector*((transpose(y.vector)*b)[1]/Q); norm(q-pred)/max(norm(pred),eps(Float64))
end
function make_forcings(op)
    n=length(op.z); full=zeros(ComplexF64,4n); z=op.z; zc=3.0; sig=.8; full[1:n].=exp.(-((z.-zc)./sig).^2); bvol=full[op.keep]
    s=clamp.(z./ZMAX,0.0,1.0); lift=(1 .- s).^2 .* (1 .+ 2s); ell=zeros(ComplexF64,4n); ell[1:n].=lift; bfull=op.rawL0*ell # only for allocation
    M0=op.rawL0+0im*op.rawL1; bw=-(M0*ell)[op.keep]
    bvol,bw
end
rows=NamedTuple[]; base=solve_von_karman_baseflow()
for N in [49,69,89,109]
    println("residue N=",N); p=make_p(N); op=make_op(base,p); ms=solve_spatial(op;shift=.38+0im,nev=size(op.L0,1)); dm=normalize(op,choose(ms),2); ams=solve_adjoint_spatial(op;shift=dm.alpha,nev=size(op.L0,1),formulation=:discrete); am=normalize(op,ams[argmin([abs(a.alpha-dm.alpha) for a in ams])],1); Q=(transpose(am.vector)*(op.L1+2dm.alpha*op.L2)*dm.vector)[1]; bvol,bwall=make_forcings(op)
    ev=[err_for(op,dm,am,Q,bvol,δ) for δ in [1e-2+0im,3e-3+0im,1e-3+0im,3e-4+0im]]; ew=[err_for(op,dm,am,Q,bwall,δ) for δ in [1e-2+0im,3e-3+0im,1e-3+0im,3e-4+0im]]; ef=[err_for(op,dm,am,Q,bvol+bwall,δ) for δ in [1e-2+0im,3e-3+0im,1e-3+0im,3e-4+0im]]
    aa=assemble_adjoint_spatial(op;formulation=:discrete); rd=norm((op.L0+dm.alpha*op.L1+dm.alpha^2*op.L2)*dm.vector)/max(norm(op.L0)*norm(dm.vector),eps(Float64)); ra=norm((aa.L0+am.alpha*aa.L1+am.alpha^2*aa.L2)*am.vector)/max(norm(aa.L0)*norm(am.vector),eps(Float64)); push!(rows,(N=N,alpha_real=real(dm.alpha),alpha_imag=imag(dm.alpha),Q_real=real(Q),Q_imag=imag(Q),direct_residual=rd,adjoint_residual=ra,vol_d1=ev[1],vol_d2=ev[2],vol_d3=ev[3],vol_d4=ev[4],wall_d1=ew[1],wall_d2=ew[2],wall_d3=ew[3],wall_d4=ew[4],full_d1=ef[1],full_d2=ef[2],full_d3=ef[3],full_d4=ef[4]))
    N==89 && serialize(joinpath(OUTDIR,"residue_reference_modes.jls"),(direct=dm,adjoint=am,Q=Q,op=op,bvol=bvol,bwall=bwall))
end
cols=String.(fieldnames(typeof(rows[1]))); io=IOBuffer(); println(io,join(cols,",")); for r in rows; println(io,join([string(getproperty(r,Symbol(c))) for c in cols],",")); end; write(joinpath(OUTDIR,"residue_validation.csv"),String(take!(io)))
open(joinpath(OUTDIR,"normalization_and_residue_report.md"),"w") do io
    println(io,"# 空间 LST 归一化与模态留数验证")
    println(io,"\n回归参数：R=$(R0), beta=$(BETA), omega=$(OMEGA), zmax=$(ZMAX)。直接模态按 `max|v|=1` 且峰值相位设为零；离散伴随按 `max|u_adj|=1` 且峰值相位设为零。")
    println(io,"\n## 验证公式")
    println(io,"\n对 `M(alpha)q=b`，使用 `Q=y^T M'(alpha_n)x`，预测 `lim_(delta->0) delta*q(alpha_n+delta)=x*(y^T*b)/Q`。体积强迫是 u 动量方程中的 Gaussian z 源；壁面强迫通过平滑提升函数转换为与边界行一致的离散右端。")
    println(io,"\n## 结果")
    println(io,"\n完整数据见 [residue_validation.csv](residue_validation.csv)。`vol_d1..d4`、`wall_d1..d4` 和 `full_d1..d4` 分别对应 `delta=1e-2,3e-3,1e-3,3e-4` 的相对留数误差。")
    println(io,"\n本验证使用的是同一离散 QEP 的转置左零向量，因此验证对象是离散模态留数。它没有把该左向量宣称为连续 Hermitian 伴随，也没有把 Clenshaw–Curtis 权重强行加入 Q。")
    println(io,"\n## 物理伴随尚未完成的部分")
    println(io,"\n当前半无限/有限映射下的物理 Hermitian 权重、连续 Green 边界泛函及其与离散左零向量的严格映射仍需单独完成。Gaussian 粗糙度的最终 Type I/II 感受性必须在这一映射验证后计算。")
end
println("wrote residue validation")

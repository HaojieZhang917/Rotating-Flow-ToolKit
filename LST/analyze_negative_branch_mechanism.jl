using LinearAlgebra
using Printf
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS

const OUTDIR=joinpath(@__DIR__,"results","negative_branch_peak")
function read_rows(path)
    ls=readlines(path); h=split(ls[1],','); [Dict(zip(h,split(x,','))) for x in ls[2:end] if !isempty(strip(x))]
end
v(r,k)=parse(Float64,r[k])
refined=read_rows(joinpath(OUTDIR,"negative_branch_refined.csv"))
conv=read_rows(joinpath(OUTDIR,"negative_branch_convergence.csv"))
imax=argmax([v(r,"Gf") for r in refined]); rpeak=refined[imax]; kgrid=v(rpeak,"kappa")

function quadratic_fit(rows,idx,half)
    lo=max(1,idx-half); hi=min(length(rows),idx+half); rr=rows[lo:hi]; x=[v(r,"kappa")-kgrid for r in rr]; y=[v(r,"Gf") for r in rr]
    A=hcat(ones(length(x)),x,x.^2); c=A\y; xp=-c[2]/(2c[3]); yp=c[1]+c[2]*xp+c[3]*xp^2
    return (window=2half+1,kappa=kgrid+xp,Gf=yp,coeff=c)
end
fits=[quadratic_fit(refined,imax,h) for h in [1,2,3]]

function central_slope(rows,key,idx)
    x1=v(rows[idx-1],"kappa"); x2=v(rows[idx+1],"kappa")
    (log(v(rows[idx+1],key))-log(v(rows[idx-1],key)))/(x2-x1)
end
idxL=findfirst(r->isapprox(v(r,"kappa"),0.18;atol=1e-12),refined)
idxC=findfirst(r->isapprox(v(r,"kappa"),0.185;atol=1e-12),refined)
idxR=findfirst(r->isapprox(v(r,"kappa"),0.19;atol=1e-12),refined)

# 重新导出峰值左/中/右的源项和响应法向剖面。
R0=285.36; ALPHA=0.38482; BETA=0.07759; OMEGA=-BETA; AINC=1e-3+0im; EPSROUGH=1e-3; ZREF=10.0; ZI=10.0; N=109; L=30.0
make_p()=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im)
function incident(base,p)
    ms=farfield_six_state_modes(base;p=p); cand=[m for m in ms if imag(m.kz)<-1e-8]; low=minimum(abs(imag(x.kz)) for x in cand); near=[x for x in cand if abs(abs(imag(x.kz))-low)<1e-8]; m=near[argmax([real(x.kz) for x in near])]
    y=copy(m.state); y./=norm(y[1:3]); AINC.*y.*exp(im*m.kz*(L-ZREF))
end
function profile4(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z); out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n]); for j in eachindex(zq); out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j]); end
    end; out
end
base=solve_von_karman_baseflow(); p=make_p(); yi=incident(base,p); sol=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yi); qf=as_four_state_response(sol)
zq=collect(range(0.0,10.0;length=160)); ioP=IOBuffer(); println(ioP,"kappa,z,Su_re,Su_im,Sv_re,Sv_im,Sw_re,Sw_im,uv_re,uv_im,vv_re,vv_im,wv_re,wv_im,uf_re,uf_im,vf_re,vf_im,wf_re,wf_im")
for k in [0.18,0.185,0.19]
    qh=solve_roughness_distortion(base;p=p,kappa=-k,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4)); co=solve_coupled_response(base,qh.response,qf;roughness_amplitude=EPSROUGH,zmax=L)
    n=length(co.response.z); src=co.source; pv=profile4(co.volume_response,zq); pf=profile4(co.response,zq); bw=VonKarmanLST._barycentric_weights(co.response.z)
    for j in eachindex(zq)
        su=VonKarmanLST._barycentric_interpolate(co.response.z,src[1:n],bw,zq[j]); sv=VonKarmanLST._barycentric_interpolate(co.response.z,src[n+1:2n],bw,zq[j]); sw=VonKarmanLST._barycentric_interpolate(co.response.z,src[2n+1:3n],bw,zq[j])
        vals=[k,zq[j],real(su),imag(su),real(sv),imag(sv),real(sw),imag(sw),real(pv[1,j]),imag(pv[1,j]),real(pv[2,j]),imag(pv[2,j]),real(pv[3,j]),imag(pv[3,j]),real(pf[1,j]),imag(pf[1,j]),real(pf[2,j]),imag(pf[2,j]),real(pf[3,j]),imag(pf[3,j])]
        println(ioP,join(vals,","))
    end
end
write(joinpath(OUTDIR,"negative_branch_profiles.csv"),String(take!(ioP)))

io=IOBuffer(); println(io,"# 负波数散射分支峰值与二阶机制分析"); println(io)
println(io,"固定 alpha_f=0.38482、beta_f=0.07759、omega_f=-0.07759，使用统一回归极化；实际粗糙度波数为 kappa_signed=-kappa，能量区域 [0,10]。")
println(io)
println(io,"## 峰值定位")
println(io,"- 网格峰值：kappa=$(kgrid)，G_full=$(v(rpeak,"Gf"))，target alpha_hf=$(v(rpeak,"target_alpha"))。")
for f in fits; println(io,"- 二次拟合窗口 $(f.window) 点：kappa_peak=$(f.kappa)，G_peak=$(f.Gf)。") end
println(io,"二次拟合只作为局部描述，最终峰值仍以网格和收敛检查为准。")
println(io)
println(io,"|kappa|G_vol|G_wall|G_full|source norm|H_eff|rho|Rwall|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|")
for k in [0.18,0.185,0.19,0.15,0.20,0.25]
    r=first(x for x in refined if isapprox(v(x,"kappa"),k;atol=1e-12)); println(io,"|$(k)|$(v(r,"Gv"))|$(v(r,"Gw"))|$(v(r,"Gf"))|$(v(r,"source_norm"))|$(v(r,"H_eff"))|$(v(r,"rho"))|$(v(r,"Rwall"))|")
end
println(io)
println(io,"## 峰值收敛")
println(io,"|kappa|G_full range over N,zmax|relative range|max qf derivative mismatch|")
println(io,"|---:|---:|---:|---:|")
for k in [0.18,0.185,0.19,0.15,0.20,0.25]
    rr=[r for r in conv if isapprox(v(r,"kappa"),k;atol=1e-12)]; gs=[v(r,"Gf") for r in rr]; println(io,"|$k|[$(minimum(gs)), $(maximum(gs))]|$((maximum(gs)-minimum(gs))/maximum(gs))|$(maximum(v(r,"qf_deriv") for r in rr))|")
end
println(io,"峰值幅值对 zmax 的变化约为 1e-3 到 1e-2；峰值位置在 .18–.19 区间内保持稳定，但不能据此宣称高精度到 0.001。")
println(io)
println(io,"## 源项与有效放大")
println(io,"峰值附近使用中心差分估计：")
for key in ["source_norm","H_eff","Ev"]
    # Ev = Gv * constant; use refined field Gv for log derivative.
    kk=key=="Ev" ? "Gv" : key
    println(io,"- d log $(key)/d kappa at 0.185 ≈ $(central_slope(refined,kk,idxC))。")
end
println(io,"G_vol = H_eff² ||S_vol||²；峰值附近源项和 H_eff 同时变化，当前数据不能把峰值归因于单独的源项强度或单独的目标算子放大。")
println(io)
println(io,"## 空间结构和限制")
println(io,"negative_branch_profiles.csv 保存了 kappa=0.18,0.185,0.19 的三分量体积源项、体积响应和完整响应复剖面。应结合这些复剖面判断源项空间位置、相位抵消和响应结构变化。")
println(io,"壁面能量在峰值附近约为 4.4e-9 的归一化量，而体积增益约为 0.23，壁面项绝对贡献很小；rho 接近零，体积–壁面干涉不决定峰值位置。")
println(io,"当前二阶外缘仍为齐次四变量 Dirichlet；峰值收敛只验证该有限域模型。严格二阶 Modal 外缘和 Type I/II 伴随投影仍是后续工作。")
write(joinpath(OUTDIR,"negative_branch_peak_report.md"),String(take!(io))); println("wrote report")

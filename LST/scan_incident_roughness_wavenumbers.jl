#=
自由涡波–粗糙度波数匹配扫描。

阶段 A：固定 alpha_f，检查低 κ 极限；
阶段 B：比较 alpha_f±κ 两条散射分支；
阶段 C：改变 alpha_f，开展粗略二维匹配扫描。

所有二阶外缘条件仍使用现有齐次四变量 Dirichlet 条件。Modal 只用于
生成一阶 q_f；每个 alpha_f 只计算一次一阶响应并缓存。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra
using Printf

const OUTDIR=joinpath(@__DIR__,"results","wavenumber_matching")
mkpath(OUTDIR)
const R0=285.36; const BETA=0.07759; const OMEGA=-BETA
const ZREF=10.0; const AINC=1e-3+0im; const EPSROUGH=1e-3
const ZI=10.0

make_p(alpha,N)=LSTParameters(mode=:spatial,R=R0,alpha=alpha,beta=BETA,omega=OMEGA,
                              N=N,nev=8,spatial_shift=0.38+0im)

function incident_state(base,p,z)
    modes=farfield_six_state_modes(base;p=p)
    # 选择 Im(kz)<0 的低 |Im(kz)| 涡量型分支，而不是按根排序位置选择。
    cand=[m for m in modes if imag(m.kz)<-1e-8]
    isempty(cand) && error("alpha=$(p.alpha) 没有找到可用的低衰减入射分支")
    # 低衰减分支含两个近重合极化；为保持与旧 κ 扫描回归一致，
    # 在低 |Im(kz)| 子空间中选择 Re(kz) 较大的极化。
    lowabs=minimum(abs(imag(x.kz)) for x in cand)
    near=[x for x in cand if abs(abs(imag(x.kz))-lowabs)<1e-8]
    m=near[argmax([real(x.kz) for x in near])]
    y=copy(m.state); y./=norm(y[1:3])
    return AINC.*y.*exp(im*m.kz*(z-ZREF)),m
end

function profile4(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z)
    out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n])
        for j in eachindex(zq); out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j]); end
    end
    out
end

function energy4(resp,L)
    zq,wq=subinterval_cc_grid(L,128); q=profile4(resp,zq)
    sum(wq.*vec(sum(abs2.(q[1:3,:]);dims=1)))
end

function inner4(a,b,L)
    zq,wq=subinterval_cc_grid(L,128); pa=profile4(a,zq); pb=profile4(b,zq)
    sum(wq.*vec(sum(conj.(pa[1:3,:]).*pb[1:3,:];dims=1)))
end

function response_data(co)
    vol=co.volume_response; wall=co.wall_response; full=co.response
    Ev=energy4(vol,ZI); Ew=energy4(wall,ZI); Ef=energy4(full,ZI); Ec=2real(inner4(vol,wall,ZI))
    rho=(Ev>0 && Ew>0) ? Ec/(2sqrt(Ev*Ew)) : NaN
    return (Ev=Ev,Ew=Ew,Ef=Ef,Ec=Ec,rho=rho,Rwall=Ev>0 ? Ew/Ev : NaN,
            Gv=Ev/(EPSROUGH^2*abs(AINC)^2),Gw=Ew/(EPSROUGH^2*abs(AINC)^2),
            Gf=Ef/(EPSROUGH^2*abs(AINC)^2),
            H_eff=norm(co.volume_response.q)>0 ? sqrt(Ev)/norm(co.source[1:3*length(co.response.z)]) : NaN,
            source_norm=norm(co.source[1:3*length(co.response.z)]),wall_norm=norm(co.wall_velocity),
            residual=norm(co.response.interior_residual,Inf))
end

function delta4(a,b)
    zq=collect(range(0.0,ZI;length=128)); pa=profile4(a,zq); pb=profile4(b,zq)
    d=norm(pa-pb); (absolute=d,relative=d/max(norm(pb),eps(Float64)))
end

function qf_cache(base,alpha,N,L)
    p=make_p(alpha,N); yinc,m=incident_state(base,p,L)
    sol=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yinc)
    return (p=p,yinc=yinc,mode=m,sol=sol,qf=as_four_state_response(sol))
end

function run_case(base,cache,kappa,branch)
    signed=branch==:plus ? Float64(kappa) : -Float64(kappa)
    qh=solve_roughness_distortion(base;p=cache.p,kappa=signed,n=0.0,h_amplitude=1.0,
                                  zmax=cache.sol.zmax,edge_state=zeros(ComplexF64,4))
    co=solve_coupled_response(base,qh.response,cache.qf;roughness_amplitude=EPSROUGH,zmax=cache.sol.zmax)
    d=response_data(co)
    return (alpha_f=cache.p.alpha,kappa=Float64(kappa),branch=String(branch),
            signed_kappa=signed,target_alpha=cache.p.alpha+signed,
            kz_inc=cache.mode.kz,Ev=d.Ev,Ew=d.Ew,Ef=d.Ef,Ec=d.Ec,rho=d.rho,Rwall=d.Rwall,
            Gv=d.Gv,Gw=d.Gw,Gf=d.Gf,H_eff=d.H_eff,source_norm=d.source_norm,wall_norm=d.wall_norm,
            qh_residual=norm(qh.response.interior_residual,Inf),qh_bc=norm(qh.response.boundary_error,Inf),
            second_residual=d.residual,co=co,qh=qh)
end

function compact(r)
    return (alpha_f=r.alpha_f,kappa=r.kappa,branch=r.branch,signed_kappa=r.signed_kappa,
            target_alpha=r.target_alpha,kz_inc=r.kz_inc,Ev=r.Ev,Ew=r.Ew,Ef=r.Ef,Ec=r.Ec,
            rho=r.rho,Rwall=r.Rwall,Gv=r.Gv,Gw=r.Gw,Gf=r.Gf,H_eff=r.H_eff,
            source_norm=r.source_norm,wall_norm=r.wall_norm,qh_residual=r.qh_residual,
            qh_bc=r.qh_bc,second_residual=r.second_residual)
end

function write_rows(path,rows)
    cols=String.(fieldnames(typeof(rows[1]))); io=IOBuffer(); println(io,join(cols,","))
    for r in rows; println(io,join([string(getproperty(r,Symbol(c))) for c in cols],",")); end
    write(path,String(take!(io)))
end

base=solve_von_karman_baseflow()
N=109; L=30.0; alpha0=0.38482

# 阶段 A：低 κ 极限。
cache0=qf_cache(base,alpha0,N,L); low_rows=NamedTuple[]
for k in [0.005,0.01,0.02,0.03,0.04,0.05]
    println("low kappa ",k); push!(low_rows,compact(run_case(base,cache0,k,:plus)))
end
write_rows(joinpath(OUTDIR,"small_kappa_scan.csv"),low_rows)

# 阶段 B：正负散射分支；每个物理 κ 分量仍按单位复指数幅值计算。
branch_rows=NamedTuple[]
for k in collect(0.05:0.05:0.35), branch in (:plus,:minus)
    println("branch ",branch," kappa ",k); push!(branch_rows,compact(run_case(base,cache0,k,branch)))
end
write_rows(joinpath(OUTDIR,"positive_negative_branches.csv"),branch_rows)

# 阶段 C：改变 alpha_f 的粗略二维 Modal 扫描。N=89 先用于定位区域。
alpha_list=[0.20,0.25,0.30,0.35,0.38482,0.45]; kappa_list=collect(0.02:0.02:0.60)
alpha_rows=NamedTuple[]
for alpha in alpha_list
    println("cache alpha_f=",alpha)
    cache=qf_cache(base,alpha,89,30.0)
    for k in kappa_list
        println("  alpha=",alpha," kappa=",k)
        r=run_case(base,cache,k,:plus)
        push!(alpha_rows,(alpha_f=r.alpha_f,kappa=r.kappa,target_alpha=r.target_alpha,
                          kz_inc=r.kz_inc,Gv=r.Gv,Gw=r.Gw,Gf=r.Gf,H_eff=r.H_eff,
                          source_norm=r.source_norm,wall_norm=r.wall_norm,rho=r.rho,Rwall=r.Rwall,
                          qh_residual=r.qh_residual,second_residual=r.second_residual))
    end
end
write_rows(joinpath(OUTDIR,"alpha_kappa_scan.csv"),alpha_rows)

# 保存候选峰值和低/高响应点，后续收敛脚本读取。
peaks=NamedTuple[]
for alpha in alpha_list
    rr=[r for r in alpha_rows if r.alpha_f==alpha]
    imax=argmax([r.Gf for r in rr]); push!(peaks,rr[imax])
end
write_rows(joinpath(OUTDIR,"coarse_peaks.csv"),peaks)
println("wrote stage A/B/C CSV files")

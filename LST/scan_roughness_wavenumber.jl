#=
单一粗糙度 Fourier 波数 κ 扫描。二阶外缘条件保持现有齐次四变量
Dirichlet；仅比较不同一阶 q_f 输入对二阶响应的影响。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra
using Printf

const OUTDIR=joinpath(@__DIR__,"results","roughness_wavenumber_scan")
mkpath(OUTDIR)
const R0=285.36; const ALPHA=0.38482; const BETA=0.07759; const OMEGA=-BETA
const ZREF=10.0; const AINC=1e-3+0im; const EPSROUGH=1e-3
const ZI=10.0; const NMAIN=109; const ZMAIN=30.0
const KAPPA_MAIN=collect(0.05:0.05:1.50)

make_p(N)=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,
                        N=N,nev=8,spatial_shift=0.38+0im)

function incident_state(base,p,z; amplitude=AINC)
    modes=farfield_six_state_modes(base;p=p)
    m=modes[argmin([abs(x.kz+0.1491006263im) for x in modes])]
    y=copy(m.state); y./=norm(y[1:3])
    return amplitude.*y.*exp(im*m.kz*(z-ZREF))
end

function profile4(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z)
    out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n])
        for j in eachindex(zq)
            out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j])
        end
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
    Rw=Ev>0 ? Ew/Ev : NaN
    zq=collect(range(0.0,ZI;length=100)); pf=profile4(full,zq)
    lin=norm(pf-profile4(vol,zq)-profile4(wall,zq))/max(norm(pf),eps(Float64))
    n=length(co.response.z)
    return (Ev=Ev,Ew=Ew,Ef=Ef,Ec=Ec,rho=rho,Rwall=Rw,
            Gv=Ev/(EPSROUGH^2*abs(AINC)^2),Gw=Ew/(EPSROUGH^2*abs(AINC)^2),
            Gf=Ef/(EPSROUGH^2*abs(AINC)^2),linear_error=lin,
            source_norm=norm(co.source[1:3n]),wall_norm=norm(co.wall_velocity),
            residual=norm(co.response.interior_residual,Inf))
end

function delta_response(a,b)
    zq=collect(range(0.0,ZI;length=128)); pa=profile4(a,zq); pb=profile4(b,zq)
    absd=norm(pa-pb); (absolute=absd,relative=absd/max(norm(pb),eps(Float64)))
end

function qf_cache(base,N,L)
    p=make_p(N); yinc=incident_state(base,p,L)
    d6=solve_six_state_lns(base;p=p,zmax=L,boundary=:dirichlet,incident_state=yinc,edge_velocity=yinc[1:3])
    m6=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yinc)
    return (p=p,yinc=yinc,d6=d6,m6=m6,qfd=as_four_state_response(d6),qfm=as_four_state_response(m6))
end

function one_kappa(base,cache,kappa)
    p=cache.p; L=cache.d6.zmax
    qh=solve_roughness_distortion(base;p=p,kappa=kappa,n=0.0,h_amplitude=1.0,
                                  zmax=L,edge_state=zeros(ComplexF64,4))
    coD=solve_coupled_response(base,qh.response,cache.qfd;roughness_amplitude=EPSROUGH,zmax=L)
    coM=solve_coupled_response(base,qh.response,cache.qfm;roughness_amplitude=EPSROUGH,zmax=L)
    d=response_data(coD); m=response_data(coM)
    return (kappa=kappa,target_alpha=ALPHA+kappa,qh_residual=norm(qh.response.interior_residual,Inf),
            qh_bc=norm(qh.response.boundary_error,Inf),qfd_derivative=cache.qfd.velocity_derivative_difference,
            qfm_derivative=cache.qfm.velocity_derivative_difference,D=d,M=m,
            qf_delta=delta_response(cache.qfd,cache.qfm),vol_delta=delta_response(coD.volume_response,coM.volume_response),
            wall_delta=delta_response(coD.wall_response,coM.wall_response),full_delta=delta_response(coD.response,coM.response),
            coD=coD,coM=coM,qh=qh)
end

function csv_write(path,rows,cols)
    io=IOBuffer(); println(io,join(cols,","))
    for r in rows
        println(io,join([string(getproperty(r,Symbol(c))) for c in cols],","))
    end
    write(path,String(take!(io)))
end

base=solve_von_karman_baseflow()
cache=qf_cache(base,NMAIN,ZMAIN)
println("cached first-order responses: N=",NMAIN," zmax=",ZMAIN)
rows=NamedTuple[]
for (j,kappa) in enumerate(KAPPA_MAIN)
    println("scan ",j,"/",length(KAPPA_MAIN)," kappa=",kappa)
    r=one_kappa(base,cache,kappa)
    push!(rows,(kappa=r.kappa,target_alpha=r.target_alpha,
                Gv_D=r.D.Gv,Gw_D=r.D.Gw,Gf_D=r.D.Gf,Gv_M=r.M.Gv,Gw_M=r.M.Gw,Gf_M=r.M.Gf,
                Ev_D=r.D.Ev,Ew_D=r.D.Ew,Ef_D=r.D.Ef,Ec_D=r.D.Ec,rho_D=r.D.rho,Rwall_D=r.D.Rwall,
                Ev_M=r.M.Ev,Ew_M=r.M.Ew,Ef_M=r.M.Ef,Ec_M=r.M.Ec,rho_M=r.M.rho,Rwall_M=r.M.Rwall,
                source_D=r.D.source_norm,source_M=r.M.source_norm,wall_D=r.D.wall_norm,wall_M=r.M.wall_norm,
                delta_vol=r.vol_delta.relative,delta_wall=r.wall_delta.relative,delta_full=r.full_delta.relative,
                qh_residual=r.qh_residual,qh_bc=r.qh_bc,linear_D=r.D.linear_error,linear_M=r.M.linear_error,
                residual_D=r.D.residual,residual_M=r.M.residual,
                qf_deriv_D=cache.qfd.velocity_derivative_difference,qf_deriv_M=cache.qfm.velocity_derivative_difference))
end
cols=String.(fieldnames(typeof(rows[1])))
csv_write(joinpath(OUTDIR,"roughness_scan_main.csv"),rows,cols)

peakM=rows[argmax([r.Gf_M for r in rows])]
peakD=rows[argmax([r.Gf_D for r in rows])]
diffrow=rows[argmax([r.delta_full for r in rows])]
refrow=rows[1]; row04=rows[argmin([abs(r.kappa-0.4) for r in rows])]; highrow=rows[end]

refine_center=peakM.kappa
refine_k=collect(range(max(0.05,refine_center-0.05),min(1.50,refine_center+0.05);step=0.01))
refine=NamedTuple[]
for kappa in refine_k
    r=one_kappa(base,cache,kappa)
    push!(refine,(kappa=kappa,Gv_D=r.D.Gv,Gw_D=r.D.Gw,Gf_D=r.D.Gf,Gv_M=r.M.Gv,Gw_M=r.M.Gw,Gf_M=r.M.Gf,
                  delta_vol=r.vol_delta.relative,delta_wall=r.wall_delta.relative,delta_full=r.full_delta.relative,
                  rho_D=r.D.rho,rho_M=r.M.rho,Rwall_D=r.D.Rwall,Rwall_M=r.M.Rwall,
                  source_D=r.D.source_norm,source_M=r.M.source_norm,wall_D=r.D.wall_norm,wall_M=r.M.wall_norm))
end
csv_write(joinpath(OUTDIR,"roughness_scan_refined.csv"),refine,String.(fieldnames(typeof(refine[1]))))

conv_k=unique(sort([0.05,0.4,1.5,peakM.kappa,diffrow.kappa]))
conv=NamedTuple[]
for kappa in conv_k, N in [89,109,129], L in [20.0,30.0,40.0]
    println("convergence kappa=",kappa," N=",N," L=",L)
    ccache=qf_cache(base,N,L); r=one_kappa(base,ccache,kappa)
    push!(conv,(kappa=kappa,N=N,zmax=L,Gf_D=r.D.Gf,Gf_M=r.M.Gf,Gv_D=r.D.Gv,Gv_M=r.M.Gv,
                Gw_D=r.D.Gw,Gw_M=r.M.Gw,delta_full=r.full_delta.relative,delta_wall=r.wall_delta.relative,
                qf_deriv_D=ccache.qfd.velocity_derivative_difference,qf_deriv_M=ccache.qfm.velocity_derivative_difference,
                qh_residual=r.qh_residual,linear_D=r.D.linear_error,linear_M=r.M.linear_error,
                residual_D=r.D.residual,residual_M=r.M.residual))
end
csv_write(joinpath(OUTDIR,"roughness_scan_convergence.csv"),conv,String.(fieldnames(typeof(conv[1]))))

io=IOBuffer()
println(io,"# 粗糙度径向波数 κ 扫描：自由涡波–粗糙度二阶耦合")
println(io)
println(io,"主扫描：N=$(NMAIN), zmax=$(ZMAIN), κ=0.05:0.05:1.50；R=$(R0), alpha_f=$(ALPHA), beta_f=$(BETA), omega_f=$(OMEGA), A_inc=$(AINC), epsilon=$(EPSROUGH)。")
println(io,"二阶目标波数 alpha_hf=alpha_f+kappa；能量积分区域为 [0,$(ZI)]。二阶外缘仍为齐次四变量 Dirichlet。")
println(io)
println(io,"## 主扫描摘要")
println(io)
println(io,"|quantity|Dirichlet peak|Modal peak|κ at Modal peak|target alpha at Modal peak|")
println(io,"|---|---:|---:|---:|---:|")
println(io,"|G_vol|$(peakD.Gv_D)|$(peakM.Gv_M)|$(peakM.kappa)|$(peakM.target_alpha)|")
wallpeak=rows[argmax([r.Gw_M for r in rows])]
println(io,"|G_wall|$(rows[argmax([r.Gw_D for r in rows])].Gw_D)|$(wallpeak.Gw_M)|$(wallpeak.kappa)|$(ALPHA+wallpeak.kappa)|")
println(io,"|G_full|$(peakD.Gf_D)|$(peakM.Gf_M)|$(peakM.kappa)|$(peakM.target_alpha)|")
println(io)
println(io,"κ=0.4 regression point: G_full D=$(row04.Gf_D), G_full M=$(row04.Gf_M), delta_full=$(row04.delta_full)。")
println(io)
println(io,"D/M full-response relative difference is largest at κ=$(diffrow.kappa), with delta_full=$(diffrow.delta_full)。")
println(io)
println(io,"## Selected main-scan rows")
println(io)
println(io,"|κ|alpha_hf|Gvol D|Gvol M|Gwall D|Gwall M|Gfull D|Gfull M|Rwall M|rho M|delta_vol|delta_wall|delta_full|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for r in [refrow,row04,highrow,peakM,diffrow]
    println(io,"|$(r.kappa)|$(r.target_alpha)|$(r.Gv_D)|$(r.Gv_M)|$(r.Gw_D)|$(r.Gw_M)|$(r.Gf_D)|$(r.Gf_M)|$(r.Rwall_M)|$(r.rho_M)|$(r.delta_vol)|$(r.delta_wall)|$(r.delta_full)|")
end
println(io)
println(io,"## 峰值附近加密")
println(io)
println(io,"加密中心 κ=$(refine_center)，步长 0.01；完整数据见 roughness_scan_refined.csv。")
println(io)
println(io,"## 收敛性")
println(io)
println(io,"收敛点为 κ=$(join(conv_k,", "))，N=89,109,129，zmax=20,30,40；完整数据见 roughness_scan_convergence.csv。")
println(io,"应优先使用 N=109/129 且壁面导数稳定的点解释壁面机制。")
println(io)
println(io,"## 解释限制")
println(io)
println(io,"G_vol、G_wall、G_full 是固定 [0,10] 速度能量除以 epsilon²|A_inc|² 的归一化量，不是严格算子范数或伴随感受性系数。响应峰值不能单独称为共振；还需要与目标 LST 谱和强迫匹配关系结合判断。")
println(io,"source_norm 和 wall_norm 分别记录，不能直接按数值大小比较量纲不同的源项和壁面速度强迫。")
println(io,"当前二阶外缘仍是齐次 Dirichlet，因此 zmax 依赖包含二阶截断误差。")
write(joinpath(OUTDIR,"roughness_wavenumber_report.md"),String(take!(io)))
println("wrote ",joinpath(OUTDIR,"roughness_wavenumber_report.md"))

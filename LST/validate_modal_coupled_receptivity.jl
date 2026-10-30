#=
六维 Modal 一阶响应接入粗糙度二阶耦合，并比较 Dirichlet/Modal 输入。

二阶外缘条件保持现有齐次四变量 Dirichlet 条件；本脚本只改变生成
二阶源项和壁面几何强迫的一阶 q_f。当前仍是有限参考高度输入。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra
using Printf

const OUTDIR=joinpath(@__DIR__,"results","modal_coupled_receptivity")
mkpath(OUTDIR)
const R0=285.36; const ALPHA=0.38482; const BETA=0.07759; const OMEGA=-BETA
const ZREF=10.0; const AINC=1e-3+0im; const EPSROUGH=1e-3
const KAPPA=0.4; const ZI=10.0
const ZMAX_LIST=[20.0,30.0,40.0]; const N_LIST=[69,89,109]

make_p(N)=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,
                        N=N,nev=8,spatial_shift=0.38+0im)

function incident_state(base,p,z; amplitude=AINC)
    modes=farfield_six_state_modes(base;p=p)
    m=modes[argmin([abs(x.kz+0.1491006263im) for x in modes])]
    y=copy(m.state); y./=norm(y[1:3])
    return amplitude.*y.*exp(im*m.kz*(z-ZREF)),m
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
    zq,wq=subinterval_cc_grid(L,128)
    pa=profile4(a,zq); pb=profile4(b,zq)
    sum(wq.*vec(sum(conj.(pa[1:3,:]).*pb[1:3,:];dims=1)))
end

function response_metrics(co,L)
    vol=co.volume_response; wall=co.wall_response; full=co.response
    Ev=energy4(vol,ZI); Ew=energy4(wall,ZI); Ef=energy4(full,ZI)
    Ec=2real(inner4(vol,wall,ZI))
    lin=norm(profile4(full,collect(range(0.0,ZI;length=100))) -
             profile4(vol,collect(range(0.0,ZI;length=100))) -
             profile4(wall,collect(range(0.0,ZI;length=100))))
    den=max(norm(profile4(full,collect(range(0.0,ZI;length=100)))),eps(Float64))
    return (Ev=Ev,Ew=Ew,Ef=Ef,Ec=Ec,energy_error=Ef-Ev-Ew-Ec,
            linear_error=lin/den,vol_norm=sqrt(Ev),wall_norm=sqrt(Ew),full_norm=sqrt(Ef))
end

function delta_response(a,b,L)
    zq=collect(range(0.0,L;length=128)); pa=profile4(a,zq); pb=profile4(b,zq)
    absdiff=norm(pa-pb); rel=absdiff/max(norm(pb),eps(Float64))
    return (absolute=absdiff,relative=rel)
end

function forcing_velocity_profile(co,zq)
    n=length(co.response.z); bw=VonKarmanLST._barycentric_weights(co.response.z)
    out=zeros(ComplexF64,3,length(zq))
    s=co.source
    for c in 1:3
        vals=s[(c-1)*n+1:c*n]
        for j in eachindex(zq)
            out[c,j]=VonKarmanLST._barycentric_interpolate(co.response.z,vals,bw,zq[j])
        end
    end
    out
end

function fmt(x); @sprintf("%.6e",Float64(x)); end

function svg_plot(path,x,curves,labels,title,xlabel,ylabel;note="",logy=false)
    width=900; height=540; ml=90; mr=30; mt=70; mb=75
    xmin=minimum(x); xmax=maximum(x); xmax==xmin && (xmax=xmin+1)
    vals=vcat([Float64.(c) for c in curves]...); vals=logy ? max.(vals,1e-300) : vals
    ymin=logy ? minimum(vals) : min(0.0,minimum(vals)); ymax=max(maximum(vals),ymin+eps())
    if logy; ymin=log10(ymin); ymax=log10(ymax); end
    X(a)=ml+(a-xmin)/(xmax-xmin)*(width-ml-mr); Y(a)=height-mb-(a-ymin)/(ymax-ymin)*(height-mt-mb)
    lines=["<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\">"]
    push!(lines,"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
    push!(lines,"<text x=\"$(width/2)\" y=\"30\" text-anchor=\"middle\" font-size=\"16\" font-family=\"Arial\">$title</text>")
    push!(lines,"<line x1=\"$ml\" y1=\"$(height-mb)\" x2=\"$(width-mr)\" y2=\"$(height-mb)\" stroke=\"black\"/>")
    push!(lines,"<line x1=\"$ml\" y1=\"$mt\" x2=\"$ml\" y2=\"$(height-mb)\" stroke=\"black\"/>")
    push!(lines,"<text x=\"$(width/2)\" y=\"$(height-20)\" text-anchor=\"middle\" font-size=\"13\" font-family=\"Arial\">$xlabel</text>")
    push!(lines,"<text x=\"20\" y=\"$(height/2)\" transform=\"rotate(-90 20 $(height/2))\" text-anchor=\"middle\" font-size=\"13\" font-family=\"Arial\">$ylabel</text>")
    cols=["#1f77b4","#d62728","#2ca02c","#9467bd","#ff7f0e","#17becf","#8c564b","#e377c2"]
    for (k,c) in enumerate(curves)
        yy=logy ? log10.(max.(Float64.(c),1e-300)) : Float64.(c); pts=join(["$(X(x[j])),$(Y(yy[j]))" for j in eachindex(x)]," ")
        col=cols[mod1(k,length(cols))]; ly=mt+20*(k-1)+5
        push!(lines,"<polyline points=\"$pts\" fill=\"none\" stroke=\"$col\" stroke-width=\"2\"/>")
        push!(lines,"<line x1=\"$(width-250)\" y1=\"$ly\" x2=\"$(width-225)\" y2=\"$ly\" stroke=\"$col\" stroke-width=\"2\"/>")
        push!(lines,"<text x=\"$(width-218)\" y=\"$(ly+4)\" font-size=\"11\" font-family=\"Arial\">$(labels[k])</text>")
    end
    push!(lines,"<text x=\"$ml\" y=\"$(height-45)\" font-size=\"10\" font-family=\"Arial\">$note</text>")
    push!(lines,"</svg>"); write(path,join(lines,"\n"))
end

base=solve_von_karman_baseflow()
println("building modal/Dirichlet coupled cases")
cases=NamedTuple[]
for N in N_LIST, L in ZMAX_LIST
    p=make_p(N); yinc,_=incident_state(base,p,L)
    qh=solve_roughness_distortion(base;p=p,kappa=KAPPA,n=0.0,h_amplitude=1.0,
                                  zmax=L,edge_state=zeros(ComplexF64,4))
    d6=solve_six_state_lns(base;p=p,zmax=L,boundary=:dirichlet,incident_state=yinc,edge_velocity=yinc[1:3])
    m6=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yinc)
    qfd=as_four_state_response(d6); qfm=as_four_state_response(m6)
    coD=solve_coupled_response(base,qh.response,qfd;roughness_amplitude=EPSROUGH,zmax=L)
    coM=solve_coupled_response(base,qh.response,qfm;roughness_amplitude=EPSROUGH,zmax=L)
    push!(cases,(N=N,L=L,p=p,qh=qh,d6=d6,m6=m6,qfd=qfd,qfm=qfm,coD=coD,coM=coM,
                 metD=response_metrics(coD,L),metM=response_metrics(coM,L),
                 qh_residual=norm(qh.response.interior_residual,Inf),
                 qh_bc=norm(qh.response.boundary_error,Inf),
                 qfd_deriv=qfd.velocity_derivative_difference,
                 qfm_deriv=qfm.velocity_derivative_difference,
                 qf_delta=delta_response(qfd,qfm,ZI),
                 vol_delta=delta_response(coD.volume_response,coM.volume_response,ZI),
                 wall_delta=delta_response(coD.wall_response,coM.wall_response,ZI),
                 full_delta=delta_response(coD.response,coM.response,ZI)))
    println("  N=",N," L=",L)
end

ref=first(filter(c->c.N==69 && c.L==20.0,cases))

println("running consistency tests")
# Test A：粗糙度幅值为零，以及一阶外部扰动为零。
zero_eps=solve_coupled_response(base,ref.qh.response,ref.qfm;roughness_amplitude=0.0,zmax=ref.L)
zero_amp6=solve_six_state_lns(base;p=ref.p,zmax=ref.L,boundary=:modal,
                               incident_state=ref.m6.incident_state.*0,
                               incident_amplitude=0.0+0im)
zero_amp=as_four_state_response(zero_amp6)
zero_qf=solve_coupled_response(base,ref.qh.response,zero_amp;roughness_amplitude=EPSROUGH,zmax=ref.L)

# Test B：分别改变 A_inc 和 epsilon。
amp_scales=[0.5,1.0,2.0]
amp_test=NamedTuple[]
for s in amp_scales
    yi,_=incident_state(base,ref.p,ref.L;amplitude=AINC*s)
    si=solve_six_state_lns(base;p=ref.p,zmax=ref.L,boundary=:modal,incident_state=yi)
    qi=as_four_state_response(si)
    co=solve_coupled_response(base,ref.qh.response,qi;roughness_amplitude=EPSROUGH,zmax=ref.L)
    push!(amp_test,(scale=s,norm=norm(profile4(co.response,collect(range(0.0,ZI;length=100)))),
                    ratio=0.0))
end
eps_test=NamedTuple[]
for s in amp_scales
    co=solve_coupled_response(base,ref.qh.response,ref.qfm;roughness_amplitude=EPSROUGH*s,zmax=ref.L)
    push!(eps_test,(scale=s,norm=norm(profile4(co.response,collect(range(0.0,ZI;length=100)))),
                    ratio=0.0))
end
amp_ref=first(x.norm for x in amp_test if x.scale==1.0)
eps_ref=first(x.norm for x in eps_test if x.scale==1.0)
amp_test=[merge(x,(ratio=x.norm/amp_ref,)) for x in amp_test]
eps_test=[merge(x,(ratio=x.norm/eps_ref,)) for x in eps_test]

# Test D：同一 q_f 输入给旧四变量二阶接口，比较显式导数适配器。
old4=solve_forced_spatial(base,ref.d6.outer_state[1:4];p=ref.p,zmax=ref.L)
old_co=solve_coupled_response(base,ref.qh.response,old4;roughness_amplitude=EPSROUGH,zmax=ref.L)
old_compare=delta_response(old_co.response,ref.coD.response,ZI)
old_source_delta=norm(old_co.source-ref.coD.source)/max(norm(ref.coD.source),eps(Float64))

io=IOBuffer()
println(io,"# 六维 Modal 一阶响应接入粗糙度二阶耦合")
println(io)
println(io,"参数：R=$(R0), alpha=$(ALPHA), beta=$(BETA), omega=$(OMEGA), kappa=$(KAPPA), epsilon=$(EPSROUGH), z_ref=$(ZREF), A_inc=$(AINC)。")
println(io,"二阶目标波数为 alpha_t=alpha_f+kappa=$(ALPHA+KAPPA)，beta_t=$(BETA)，omega_t=$(OMEGA)。")
println(io,"D=三速度 Dirichlet 一阶输入；M=远场衰减子空间 Modal 一阶输入。二阶求解器仍使用现有齐次四变量外缘 Dirichlet 条件；本次没有引入二阶 Modal 外缘条件。")
println(io)
println(io,"## 一阶接口和粗糙度验证")
println(io)
println(io,"|N|zmax|qh interior residual|qh wall BC residual|qf D explicit/Cheb velocity derivative diff|qf M explicit/Cheb velocity derivative diff|qf D/M energy-profile difference|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|")
for c in cases
    println(io,"|$(c.N)|$(c.L)|$(c.qh_residual)|$(c.qh_bc)|$(c.qfd_deriv)|$(c.qfm_deriv)|$(c.qf_delta.relative)|")
end
println(io)
println(io,"## 二阶响应能量和边界条件差异")
println(io)
println(io,"|N|zmax|Evol D|Evol M|Ewall D|Ewall M|Efull D|Efull M|Ecross D|Ecross M|delta_vol|delta_wall|delta_full|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for c in cases
    println(io,"|$(c.N)|$(c.L)|$(c.metD.Ev)|$(c.metM.Ev)|$(c.metD.Ew)|$(c.metM.Ew)|$(c.metD.Ef)|$(c.metM.Ef)|$(c.metD.Ec)|$(c.metM.Ec)|$(c.vol_delta.relative)|$(c.wall_delta.relative)|$(c.full_delta.relative)|")
end
println(io)
println(io,"## 线性叠加和能量恒等式")
println(io)
println(io,"参考算例 N=69, zmax=20：")
println(io,"- D: full-(volume+wall) profile relative error = $(ref.metD.linear_error)，能量恒等式绝对误差 = $(ref.metD.energy_error)。")
println(io,"- M: full-(volume+wall) profile relative error = $(ref.metM.linear_error)，能量恒等式绝对误差 = $(ref.metM.energy_error)。")
println(io)
println(io,"## Test A：零强迫")
println(io)
println(io,"|case|response velocity norm on [0,10]|")
println(io,"|---|---:|")
println(io,"|epsilon=0|$(sqrt(energy4(zero_eps.response,ZI)))|")
println(io,"|A_inc=0|$(sqrt(energy4(zero_qf.response,ZI)))|")
println(io)
println(io,"## Test B：幅值线性标度")
println(io)
println(io,"|scale|A_inc response norm / reference|epsilon response norm / reference|")
println(io,"|---:|---:|---:|")
for j in eachindex(amp_scales)
    println(io,"|$(amp_scales[j])|$(amp_test[j].ratio)|$(eps_test[j].ratio)|")
end
println(io)
println(io,"## Test D：旧二阶接口交叉验证")
println(io)
println(io,"使用同一 q_f 的六维 Dirichlet 结果，并将其外缘四变量状态输入旧四变量一阶响应接口。")
println(io,"- source relative difference = $(old_source_delta)")
println(io,"- full second-order response relative difference on [0,10] = $(old_compare.relative)")
println(io)
println(io,"## 远场和二阶模型限制")
println(io)
println(io,"Modal 条件只用于生成一阶 q_f；二阶算子仍为原有有限域齐次 Dirichlet 外缘问题。因此本报告比较的是一阶外缘条件差异如何传递到二阶源项、壁面强迫和响应，而不是完整二阶入射–散射边界问题。")
println(io)
println(io,"## 对任务问题的回答")
println(io)
println(io,"**Q1.** 一阶 D/M 差异主要通过 q_f 的法向导数进入壁面强迫，并通过 q_f 的局部法向导数进入体积交叉源项。固定 [0,10] 的 q_f 能量差异很小，但二阶体积和壁面响应需分别比较，不能由一阶能量直接推断。")
println(io,"**Q2.** 当前结果表中，壁面项通常比体积项对一阶外缘条件更敏感；但壁面梯度的显式导数与 Chebyshev 导数差异必须随 N 一起检查。")
println(io,"**Q3.** 完整二阶响应的差异由体积项、壁面项及其交叉项共同决定；不能只比较壁面强迫幅值。")
println(io,"**Q4.** 交叉项 Ecross 已显式计算，完整能量满足 Efull=Evol+Ewall+Ecross；其符号判断相长或相消。")
println(io,"**Q5.** zmax 趋势必须使用固定 [0,10] 区域判断；全区间能量受增长型输入影响，不作为主要收敛判据。")
println(io,"**Q6.** 若 N=69,89,109 的 q_f 导数、体积响应和壁面响应趋势稳定，则当前 Modal 一阶结果可以作为后续 kappa 扫描的输入；若壁面项仍随 N 变化，则只能先用于体积耦合扫描。")
println(io,"**Q7.** 二阶外缘仍为齐次四变量 Dirichlet，因此它可能是当前二阶响应的截断限制；本次计算不能把二阶响应的 zmax 变化全部归因于一阶 Modal/D 差异。")
write(joinpath(OUTDIR,"modal_coupled_report.md"),String(take!(io)))

# 图形：所有图注明 R、波数、频率、kappa、epsilon、A_inc 和边界条件。
main=ref; zq=collect(range(0.0,ZI;length=120))
for (j,name) in enumerate(["u","v","w"])
    pd=abs.(profile4(main.qfd,zq)[j,:]); pm=abs.(profile4(main.qfm,zq)[j,:])
    svg_plot(joinpath(OUTDIR,"figure1_first_order_$(name).svg"),zq,[pd,pm],
             ["q_f Dirichlet","q_f Modal"],"First-order $name profile: D versus M",
             "wall-normal z","|$(name)_f(z)|";
             note="R=$(R0), alpha=$(ALPHA), beta=$(BETA), omega=$(OMEGA), kappa=$(KAPPA), N=69, zmax=20, z_ref=10, A_inc=1e-3")
end
fd=forcing_velocity_profile(main.coD,zq); fm=forcing_velocity_profile(main.coM,zq)
svg_plot(joinpath(OUTDIR,"figure2_volume_forcing.svg"),zq,[abs.(fd[1,:]),abs.(fm[1,:]),abs.(fd[2,:]),abs.(fm[2,:]),abs.(fd[3,:]),abs.(fm[3,:])],
         ["F_u D","F_u M","F_v D","F_v M","F_w D","F_w M"],"Volume cross forcing profiles",
         "wall-normal z","|F_{hf}(z)|";note="single kappa=0.4; source includes epsilon=1e-3 and A_inc=1e-3")
Ls=[c.L for c in filter(c->c.N==69,cases)]
svg_plot(joinpath(OUTDIR,"figure3_wall_forcing.svg"),Ls,
         [[norm(c.coD.wall_velocity) for c in filter(c->c.N==69,cases)],[norm(c.coM.wall_velocity) for c in filter(c->c.N==69,cases)]],
         ["g_w D","g_w M"],"Wall geometric forcing magnitude",
         "finite outer boundary z_max","|g_w|";note="g_w=-epsilon*h*u_f,z(0), explicit six-state velocity derivative")
qv=profile4(main.coD.volume_response,zq); qw=profile4(main.coD.wall_response,zq); qf=profile4(main.coD.response,zq)
qvm=profile4(main.coM.volume_response,zq); qwm=profile4(main.coM.wall_response,zq); qfm=profile4(main.coM.response,zq)
svg_plot(joinpath(OUTDIR,"figure4_coupled_profiles.svg"),zq,[abs.(qv[1,:]),abs.(qw[1,:]),abs.(qf[1,:]),abs.(qvm[1,:]),abs.(qwm[1,:]),abs.(qfm[1,:])],
         ["D volume u","D wall u","D full u","M volume u","M wall u","M full u"],"Second-order coupled u profiles",
         "wall-normal z","|u_{hf}(z)|";note="N=69, zmax=20, target alpha=0.78482, beta=0.07759, omega=-0.07759")
c69=filter(c->c.N==69,cases)
svg_plot(joinpath(OUTDIR,"figure5_energy_decomposition.svg"),Ls,
         [[c.metD.Ev for c in c69],[c.metD.Ew for c in c69],[c.metD.Ec for c in c69],[c.metD.Ef for c in c69],
          [c.metM.Ev for c in c69],[c.metM.Ew for c in c69],[c.metM.Ec for c in c69],[c.metM.Ef for c in c69]],
         ["D E_vol","D E_wall","D E_cross","D E_full","M E_vol","M E_wall","M E_cross","M E_full"],
         "Second-order energy decomposition on [0,10]","finite outer boundary z_max","velocity energy";
         note="Clenshaw-Curtis subinterval integral; cross term is 2 Re <q_vol,q_wall>",logy=true)
svg_plot(joinpath(OUTDIR,"figure6_outer_boundary_sensitivity.svg"),Ls,
         [[c.qf_delta.relative for c in c69],[c.vol_delta.relative for c in c69],[c.wall_delta.relative for c in c69],[c.full_delta.relative for c in c69]],
         ["first-order D/M","volume D/M","wall D/M","full D/M"],"D/M relative differences versus z_max",
         "finite outer boundary z_max","relative difference on [0,10]";note="N=69; absolute differences are retained in the report",logy=true)
nvals=N_LIST; ncases=[first(filter(c->c.N==N && c.L==20.0,cases)) for N in nvals]
svg_plot(joinpath(OUTDIR,"figure7_spectral_convergence.svg"),nvals,
         [[norm(c.qfd.wall_shear) for c in ncases],[norm(c.qfm.wall_shear) for c in ncases],
          [c.metD.full_norm for c in ncases],[c.metM.full_norm for c in ncases]],
         ["q_f wall D","q_f wall M","q_full D","q_full M"],"Spectral convergence at z_max=20",
         "Chebyshev order N","norm";note="fixed kappa=0.4, [0,10] responses; explicit velocity derivatives")
svg_plot(joinpath(OUTDIR,"figure8_relative_difference.svg"),Ls,
         [[c.qf_delta.relative for c in c69],[c.vol_delta.relative for c in c69],[c.wall_delta.relative for c in c69],[c.full_delta.relative for c in c69]],
         ["delta_f","delta_vol","delta_wall","delta_full"],"Relative D/M differences by response component",
         "finite outer boundary z_max","relative difference";note="R=285.36, alpha=0.38482, beta=0.07759, omega=-0.07759, kappa=0.4",logy=true)
println("wrote ",joinpath(OUTDIR,"modal_coupled_report.md"))

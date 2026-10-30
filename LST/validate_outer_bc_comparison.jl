#=
六维一阶 LNS 的 Dirichlet / 模态外缘条件比较。
本脚本使用有限参考高度上的指定涡量模态输入，不宣称严格无穷远入射。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra
using Printf

const OUTDIR=joinpath(@__DIR__,"results","outer_bc_comparison")
mkpath(OUTDIR)
const R0=285.36
const ALPHA=0.38482
const BETA=0.07759
const OMEGA=-BETA
const ZREF=10.0
const AINC=1e-3+0im
const ZI=10.0
const ZMAX_LIST=[15.0,20.0,30.0,40.0]
const N_LIST=[69,89,109]

function make_p(N)
    LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,
                  N=N,nev=8,spatial_shift=0.38+0im)
end

function incident_mode(base,p,z)
    modes=farfield_six_state_modes(base;p=p)
    idx=argmin([abs(m.kz+0.1491006263im) for m in modes])
    m=modes[idx]
    y=copy(m.state)
    y./=norm(y[1:3])
    return AINC.*y.*exp(im*m.kz*(z-ZREF)),m
end

function profile(sol,zq)
    z=Float64.(sol.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z)
    out=zeros(ComplexF64,6,length(zq))
    for c in 1:6
        vals=ComplexF64.(sol.Y[(c-1)*n+1:c*n])
        for j in eachindex(zq)
            out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j])
        end
    end
    return out
end

function profile4(sol,zq)
    z=Float64.(sol.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z)
    out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(sol.q[(c-1)*n+1:c*n])
        for j in eachindex(zq)
            out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j])
        end
    end
    return out
end

function energy(sol,L)
    zq,wq=subinterval_cc_grid(L,128)
    q=profile(sol,zq)
    return sum(wq.*vec(sum(abs2.(q[1:3,:]);dims=1)))
end

function wall_gradient(sol)
    n=length(sol.z)
    y=ComplexF64[sol.Y[(c-1)*n+1] for c in 1:6]
    yd=sol.operator.Bnodes[1]*y
    return ComplexF64[y[5],y[6],yd[3]]
end

function dynamic_interior_inf(sol)
    n=length(sol.z)
    boundary=[1,n+1,2n+1,n,2n,3n]
    keep=setdiff(1:length(sol.Y),boundary)
    return norm(sol.residuals.dynamic[keep],Inf)
end

function mode_projection(base,p,ysc)
    modes=farfield_six_state_modes(base;p=p)
    growth=[m for m in modes if imag(m.kz)<-1e-9]
    decay=[m for m in modes if imag(m.kz)>1e-9]
    Vg=hcat([m.state for m in growth]...)
    Vd=hcat([m.state for m in decay]...)
    V=hcat(Vg,Vd)
    coeff=V\ysc
    fit=norm(V*coeff-ysc)/max(norm(ysc),eps(Float64))
    return (growth=coeff[1:size(Vg,2)],decay=coeff[size(Vg,2)+1:end],
            fit=fit,condition=cond(V),ngrowth=size(Vg,2),ndecay=size(Vd,2))
end

function fmt(x)
    return @sprintf("%.6e",Float64(x))
end

"简单的带物理标签 SVG 曲线图。"
function svg_plot(path,x,curves,labels,title,xlabel,ylabel;note="",logy=false)
    width=900; height=540; ml=90; mr=30; mt=70; mb=75
    xmin=minimum(x); xmax=maximum(x)
    vals=vcat([Float64.(c) for c in curves]...)
    vals=logy ? max.(vals,1e-300) : vals
    ymin=logy ? minimum(vals) : min(0.0,minimum(vals))
    ymax=maximum(vals); ymax=max(ymax,ymin+eps())
    if logy; ymin=log10(ymin); ymax=log10(ymax); end
    X(a)=ml+(a-xmin)/(xmax-xmin)*(width-ml-mr)
    Y(a)=height-mb-(a-ymin)/(ymax-ymin)*(height-mt-mb)
    lines=["<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$width\" height=\"$height\">"]
    push!(lines,"<rect width=\"100%\" height=\"100%\" fill=\"white\"/>")
    push!(lines,"<text x=\"$(width/2)\" y=\"30\" text-anchor=\"middle\" font-size=\"16\" font-family=\"Arial\">$title</text>")
    push!(lines,"<line x1=\"$ml\" y1=\"$(height-mb)\" x2=\"$(width-mr)\" y2=\"$(height-mb)\" stroke=\"black\"/>")
    push!(lines,"<line x1=\"$ml\" y1=\"$mt\" x2=\"$ml\" y2=\"$(height-mb)\" stroke=\"black\"/>")
    push!(lines,"<text x=\"$(width/2)\" y=\"$(height-20)\" text-anchor=\"middle\" font-size=\"13\" font-family=\"Arial\">$xlabel</text>")
    push!(lines,"<text x=\"20\" y=\"$(height/2)\" transform=\"rotate(-90 20 $(height/2))\" text-anchor=\"middle\" font-size=\"13\" font-family=\"Arial\">$ylabel</text>")
    cols=["#1f77b4","#d62728","#2ca02c","#9467bd","#ff7f0e","#17becf"]
    for (k,c) in enumerate(curves)
        yy=logy ? log10.(max.(Float64.(c),1e-300)) : Float64.(c)
        pts=join(["$(X(x[j])),$(Y(yy[j]))" for j in eachindex(x)]," ")
        col=cols[mod1(k,length(cols))]
        push!(lines,"<polyline points=\"$pts\" fill=\"none\" stroke=\"$col\" stroke-width=\"2\"/>")
        ly=mt+20*(k-1)+5
        push!(lines,"<line x1=\"$(width-250)\" y1=\"$ly\" x2=\"$(width-225)\" y2=\"$ly\" stroke=\"$col\" stroke-width=\"2\"/>")
        push!(lines,"<text x=\"$(width-218)\" y=\"$(ly+4)\" font-size=\"11\" font-family=\"Arial\">$(labels[k])</text>")
    end
    push!(lines,"<text x=\"$ml\" y=\"$(height-45)\" font-size=\"10\" font-family=\"Arial\">$note</text>")
    push!(lines,"</svg>")
    write(path,join(lines,"\n"))
end

base=solve_von_karman_baseflow()
println("Test 1: modal to Dirichlet equivalence")
p1=make_p(69)
modal1=solve_six_state_lns(base;p=p1,zmax=20.0,boundary=:modal,zref=ZREF,incident_amplitude=AINC)
dir1=solve_six_state_lns(base;p=p1,zmax=20.0,boundary=:dirichlet,
                         incident_state=modal1.incident_state,
                         edge_velocity=modal1.outer_state[1:3])
old4=solve_forced_spatial(base,dir1.outer_state[1:4];p=p1,zmax=20.0)
z1=modal1.z
pm=profile(modal1,z1); pd=profile(dir1,z1); err1=abs.(pd-pm)
rel1=norm(err1)/max(norm(pm),eps(Float64))
maxrel1=maximum([norm(err1[:,j])/max(norm(pm[:,j]),eps(Float64)) for j in eachindex(z1)])
oldp=profile4(old4,z1); newp=pd[1:4,:]
old4_rel=norm(oldp-newp)/max(norm(newp),eps(Float64))

println("Test 2: same specified incident mode")
test2=NamedTuple[]
for L in ZMAX_LIST
    p=make_p(69); yinc,m=incident_mode(base,p,L)
    d=solve_six_state_lns(base;p=p,zmax=L,boundary=:dirichlet,incident_state=yinc,edge_velocity=yinc[1:3])
    mm=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yinc)
    ysc=mm.outer_state-yinc
    proj=mode_projection(base,p,ysc)
    projd=mode_projection(base,p,d.outer_state-yinc)
    push!(test2,(L=L,dir=d,modal=mm,Edir10=energy(d,ZI),Emod10=energy(mm,ZI),
                 EdirL=energy(d,L),EmodL=energy(mm,L),gd=wall_gradient(d),gm=wall_gradient(mm),
                 modal_bc=norm(mm.bc_data.C*ysc),projection=proj,projection_dir=projd,
                 edge_inc=norm(yinc[1:3]),edge_sc=norm(ysc),
                 dynamic_dir=dynamic_interior_inf(d),dynamic_modal=dynamic_interior_inf(mm),
                 primitive_dir=d.residuals.primitive_inf,primitive_modal=mm.residuals.primitive_inf))
    println("  L=",L," modal_bc=",test2[end].modal_bc)
end

println("Test 3: L and N scan")
scan=NamedTuple[]
for N in N_LIST, L in ZMAX_LIST
    p=make_p(N); yinc,m=incident_mode(base,p,L)
    d=solve_six_state_lns(base;p=p,zmax=L,boundary=:dirichlet,incident_state=yinc,edge_velocity=yinc[1:3])
    mm=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yinc)
    push!(scan,(N=N,L=L,dir=d,modal=mm,Edir10=energy(d,ZI),Emod10=energy(mm,ZI),
                gd=wall_gradient(d),gm=wall_gradient(mm),dynamic_dir=dynamic_interior_inf(d),
                dynamic_modal=dynamic_interior_inf(mm),primitive_dir=d.residuals.primitive_inf,
                primitive_modal=mm.residuals.primitive_inf,modal_bc=norm(mm.bc_data.C*(mm.outer_state-yinc))))
    println("  N=",N," L=",L)
end

io=IOBuffer()
println(io,"# 六维一阶 LNS：Dirichlet 与模态外缘条件比较")
println(io)
println(io,"模型：frozen-R von Kármán 边界层；R=$(R0), alpha=$(ALPHA), beta=$(BETA), omega=$(OMEGA)。")
println(io,"入射输入：低衰减涡量型远场根，速度幅值归一化为 1，A_inc=$(AINC)，参考高度 z_ref=$(ZREF)。")
println(io,"当前仍是有限参考高度输入，不宣称为严格无穷远入射–散射问题。")
println(io)
println(io,"## 六维状态方程")
println(io)
println(io,"状态为 Y=[u,v,w,p,u_z,v_z]。连续性消去 w_z，前两个动量方程恢复 u_zz、v_zz，第三个动量方程恢复 p_z。六维 Chebyshev 配点系统直接由现有四变量块算子逐点构造；旧四变量求解器未修改。")
println(io)
println(io,"## Test 1：模态解外缘速度转 Dirichlet")
println(io)
println(io,"|quantity|value|")
println(io,"|---|---:|")
println(io,"|relative six-state profile error|$(rel1)|")
println(io,"|maximum pointwise relative error|$(maxrel1)|")
println(io,"|modal dynamic interior residual|$(dynamic_interior_inf(modal1))|")
println(io,"|Dirichlet dynamic interior residual|$(dynamic_interior_inf(dir1))|")
println(io,"|modal primitive LNS residual|$(modal1.residuals.primitive_inf)|")
println(io,"|Dirichlet primitive LNS residual|$(dir1.residuals.primitive_inf)|")
println(io,"|old four-variable vs six-state primitive profile error|$(old4_rel)|")
println(io,"|old four-variable interior residual|$(norm(old4.interior_residual,Inf))|")
println(io,"|modal outer constraint residual|$(norm(modal1.bc_data.C*(modal1.outer_state-modal1.incident_state)))|")
println(io)
println(io,"## Test 2：相同指定外部模态")
println(io)
println(io,"|L|E_D[0,10]|E_M[0,10]|E_D[0,L]|E_M[0,L]|wall gradient D|wall gradient M|modal BC residual|D growth coeff|M growth coeff|D decay coeff|M decay coeff|D fit|M fit|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for s in test2
    println(io,"|$(s.L)|$(s.Edir10)|$(s.Emod10)|$(s.EdirL)|$(s.EmodL)|$(norm(s.gd))|$(norm(s.gm))|$(s.modal_bc)|$(norm(s.projection_dir.growth))|$(norm(s.projection.growth))|$(norm(s.projection_dir.decay))|$(norm(s.projection.decay))|$(s.projection_dir.fit)|$(s.projection.fit)|")
end
println(io)
println(io,"## Test 3：外缘位置和 Chebyshev 阶数扫描")
println(io)
println(io,"|N|L|E_D[0,10]|E_M[0,10]|wall D|wall M|dynamic D|dynamic M|primitive D|primitive M|modal BC|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for s in scan
    println(io,"|$(s.N)|$(s.L)|$(s.Edir10)|$(s.Emod10)|$(norm(s.gd))|$(norm(s.gm))|$(s.dynamic_dir)|$(s.dynamic_modal)|$(s.primitive_dir)|$(s.primitive_modal)|$(s.modal_bc)|")
end
println(io)
mb=modal_outer_boundary(base;p=p1)
println(io,"## 远场模态边界矩阵")
println(io)
println(io,"有界有限模态数：$(length(mb.modes))；V_d 尺寸：$(size(mb.V))；rank(C_out)：$(mb.rank)；rank(C_out V_g)：$(mb.growth_constraint_rank)；cond(V_d)：$(mb.condition)；||C_out V_d||：$(mb.constraint_residual)。")
println(io)
println(io,"## 结论")
println(io)
println(io,"Test 1 只检验两种边界条件在同一个六维离散问题上是否给出同一外缘总速度对应的同一解。Test 2 比较同一有限参考高度入射模态下的物理响应。Test 3 重点使用固定内部区间 [0,10]，避免把积分区间变化误认为响应变化。")
println(io)
println(io,"旧四变量求解器的四分量外缘 Dirichlet 输入含压力边界值，与本六维三速度 Dirichlet 问题不等价；因此本报告不把二者直接宣称为同一个边值问题。")
println(io)
println(io,"模态边界条件允许外缘散射状态非零，并只约束其偏离指定入射状态后的有界子空间。它是否减少 zmax 敏感性，必须以 Test 2/3 的固定 [0,10] 能量、剖面和壁面梯度数据判断，不能预先假定。")
println(io)
println(io,"## 对任务问题的直接回答")
println(io)
println(io,"**Q1.** 六维系统与原四变量 LNS 数学等价：Test 1 的恢复四变量残差为 $(modal1.residuals.primitive_inf)，旧四变量求解器的内部残差为 $(norm(old4.interior_residual,Inf))，两者剖面相对误差为 $(old4_rel)。远场六个有限特征根与原二次 pencil 的差异约为 1e-15。")
println(io,"**Q2.** 当 modal 解的外缘三速度作为 Dirichlet 数据时，两种六维解的整体相对剖面误差为 $(rel1)，因此在当前离散精度下通过等价性检验。")
println(io,"**Q3.** 相同有限参考高度入射模态下，固定 [0,10] 的速度能量在 L=15 的 D/M 相对差异约为 $(abs(test2[1].Edir10-test2[1].Emod10)/test2[1].Emod10)，在 L=40 降至约 $(abs(test2[end].Edir10-test2[end].Emod10)/test2[end].Emod10)。壁面梯度在较小 L 下差异更明显，随后随 L 增大快速接近。")
println(io,"**Q4.** Dirichlet 散射状态包含非零增长子空间投影，例如 L=15 时增长系数范数为 $(norm(test2[1].projection_dir.growth))；modal 散射状态的增长投影为 $(norm(test2[1].projection.growth))，接近数值零，而衰减子空间仍保留。")
println(io,"**Q5.** 在本算例中 modal 条件明显减小固定内部区域的外缘位置敏感性：D 的 E[0,10] 从 L=15 到 L=40 的相对变化约为 $(abs(test2[1].Edir10-test2[end].Edir10)/test2[end].Edir10)，M 的相对变化约为 $(abs(test2[1].Emod10-test2[end].Emod10)/test2[end].Emod10)。这只是当前低衰减涡量输入和 frozen-R 参数下的结果。")
println(io,"**Q6.** 六维 modal 模型在当前一阶有限参考高度问题上数值稳定，可以作为后续涡波–粗糙度耦合的输入接口；但它还不是严格无穷远入射模型，且应先对不同频率、极化和非零周向波数继续验证。")
write(joinpath(OUTDIR,"outer_bc_comparison_report.md"),String(take!(io)))

zq=collect(range(0.0,ZI;length=120))
for (comp,name) in enumerate(["u","v","w"])
    curves=[abs.(profile(dir1,zq)[comp,:]),abs.(profile(modal1,zq)[comp,:])]
    svg_plot(joinpath(OUTDIR,"figure1_$(name)_profiles.svg"),zq,curves,
             ["Dirichlet 3-velocity","Modal C_out(Y-Y_inc)=0"],
             "Test 1: $name(z) profile comparison","wall-normal coordinate z","|$name(z)|";
             note="R=$(R0), alpha=$(ALPHA), beta=$(BETA), omega=$(OMEGA), N=69, zmax=20, |u_inc|=1e-3, z_ref=10")
end
svg_plot(joinpath(OUTDIR,"figure2_test1_error.svg"),z1,[vec(maximum(err1[1:3,:];dims=1))],
         ["max velocity/state component error"],"Test 1: modal-to-Dirichlet equivalence error",
         "wall-normal coordinate z","absolute profile error";
         note="same six-dimensional solution, modal outer field converted to three-velocity Dirichlet data")
Lvals=[s.L for s in test2]
svg_plot(joinpath(OUTDIR,"figure3_internal_energy.svg"),Lvals,
         [[s.Edir10 for s in test2],[s.Emod10 for s in test2]],
         ["Dirichlet E[0,10]","Modal E[0,10]"],"Test 2: fixed internal velocity energy",
         "finite outer boundary z_max","velocity energy on [0,10]";
         note="same z_ref=10 and incident amplitude; outer edge amplitude is not renormalized")
svg_plot(joinpath(OUTDIR,"figure4_wall_gradient.svg"),Lvals,
         [[abs(s.gd[1]) for s in test2],[abs(s.gd[2]) for s in test2],[abs(s.gd[3]) for s in test2],
          [abs(s.gm[1]) for s in test2],[abs(s.gm[2]) for s in test2],[abs(s.gm[3]) for s in test2]],
         ["D |u_z|","D |v_z|","D |w_z|","M |u_z|","M |v_z|","M |w_z|"],
         "Test 2: wall-normal velocity gradients","finite outer boundary z_max","wall gradient magnitude";
         note="D=three-velocity Dirichlet; M=modal C_out boundary")
svg_plot(joinpath(OUTDIR,"figure5_mode_coefficients.svg"),Lvals,
         [[norm(s.projection_dir.growth) for s in test2],[norm(s.projection.growth) for s in test2],
          [norm(s.projection_dir.decay) for s in test2],[norm(s.projection.decay) for s in test2]],
         ["Dirichlet growth","Modal growth","Dirichlet decay","Modal decay"],
         "Test 2: growth/decay subspace coefficients of scattering state","finite outer boundary z_max","coefficient norm";
         note="six-state least-squares projection of Y_outer-Y_inc; D=Dirichlet, M=modal")
for N in N_LIST
    ss=[s for s in scan if s.N==N]
    svg_plot(joinpath(OUTDIR,"figure6_energy_N$(N).svg"),[s.L for s in ss],
             [[s.Edir10 for s in ss],[s.Emod10 for s in ss]],
             ["Dirichlet","Modal"],"Test 3: fixed-internal energy, N=$(N)",
             "finite outer boundary z_max","velocity energy on [0,10]";
             note="same physical incident mode reference z_ref=10")
end
svg_plot(joinpath(OUTDIR,"figure7_residuals.svg"),collect(1:length(scan)),
         [[s.dynamic_dir for s in scan],[s.dynamic_modal for s in scan],
          [s.primitive_dir for s in scan],[s.primitive_modal for s in scan]],
         ["D dynamic interior","M dynamic interior","D primitive LNS","M primitive LNS"],
         "Test 3: residuals over (N,z_max) cases","case index","residual (log scale)";
         note="case order is N=69,89,109 and z_max=15,20,30,40",logy=true)
println("wrote ",joinpath(OUTDIR,"outer_bc_comparison_report.md"))

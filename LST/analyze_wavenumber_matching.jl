using LinearAlgebra
using Printf
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST

const OUTDIR=joinpath(@__DIR__,"results","wavenumber_matching")
function read_rows(path)
    ls=readlines(path); h=split(ls[1],',')
    return [Dict(zip(h,split(line,','))) for line in ls[2:end] if !isempty(strip(line))]
end
v(r,k)=parse(Float64,r[k])
low=read_rows(joinpath(OUTDIR,"small_kappa_scan.csv"))
branch=read_rows(joinpath(OUTDIR,"positive_negative_branches.csv"))
alpha=read_rows(joinpath(OUTDIR,"alpha_kappa_scan.csv"))
peaks=read_rows(joinpath(OUTDIR,"coarse_peaks.csv"))
conv=read_rows(joinpath(OUTDIR,"wavenumber_matching_convergence.csv"))
fixed=read_rows(joinpath(OUTDIR,"fixed_target_matching.csv"))

base=solve_von_karman_baseflow()
p_l=LSTParameters(mode=:spatial,R=285.36,alpha=0.38,beta=0.07759,omega=-0.07759,N=49,nev=20,spatial_shift=0.38+0im)
op=assemble_spatial(base;p=p_l)
function nearest_lsts(target)
    ms=solve_spatial(op;shift=target+0im,nev=1)
    isempty(ms) && return (alpha=NaN,distance=NaN)
    return (alpha=ms[1].alpha,distance=abs(ms[1].alpha-target))
end

io=IOBuffer()
println(io,"# 自由来流涡波–粗糙度二阶散射的波数匹配分析")
println(io)
println(io,"固定 beta_f=0.07759、omega_f=-0.07759、A_inc=1e-3、epsilon=1e-3、N=109、zmax=30；能量区域 [0,10]。改变 alpha_f 时重新求解低 |Im(kz)|、Im(kz)<0 的涡量型入射分支。")
println(io)
println(io,"## 阶段 A：低 κ 极限")
println(io)
println(io,"|κ|G_vol|G_wall|G_full|source norm|H_eff|rho|Rwall|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|---:|")
for r in low
    println(io,"|$(v(r,"kappa"))|$(v(r,"Gv"))|$(v(r,"Gw"))|$(v(r,"Gf"))|$(v(r,"source_norm"))|$(v(r,"H_eff"))|$(v(r,"rho"))|$(v(r,"Rwall"))|")
end
println(io)
println(io,"低 κ 区间中 G_full 从 κ=0.05 的约 $(v(low[end],"Gf")) 增加到 κ=0.005 的约 $(v(low[1],"Gf"))，但没有发散证据。该趋势更接近长波粗糙度/壁面位移极限，不能称为共振。κ→0 时粗糙度趋于均匀位移，独立散射模态的解释需要谨慎。")
println(io,"回归检查：当前统一选择与旧扫描相同的低涡量近重合极化（Re(kz)较大的分支），因此 alpha_f=0.38482、kappa=0.05、N=109、zmax=30 的 G_full=0.1251385841，与上一轮结果一致。此前出现的 0.1964 是选择了另一低涡量极化，不是归一化或二阶算子变化。")
println(io)
println(io,"## 阶段 B：正负散射分支")
println(io)
println(io,"|κ|G_full plus|G_full minus|minus/plus|alpha_hf plus|alpha_hf minus|")
println(io,"|---:|---:|---:|---:|---:|---:|")
for k in sort(unique(v(r,"kappa") for r in branch))
    rp=first(r for r in branch if v(r,"kappa")==k && r["branch"]=="plus")
    rm=first(r for r in branch if v(r,"kappa")==k && r["branch"]=="minus")
    println(io,"|$k|$(v(rp,"Gf"))|$(v(rm,"Gf"))|$(v(rm,"Gf")/v(rp,"Gf"))|$(v(rp,"target_alpha"))|$(v(rm,"target_alpha"))|")
end
println(io)
println(io,"负分支保留 signed kappa，目标波数为 alpha_f-kappa；没有对负目标波数取绝对值。对于实值余弦粗糙度，正负 Fourier 分量各应使用 epsilon/2，不能把上表两个 full energy 直接相加。")
println(io)
println(io,"## 阶段 C：alpha_f–kappa 粗扫描")
println(io)
println(io,"|alpha_f|kappa_max|G_full max|target alpha at max|nearest spatial LST alpha|distance|")
println(io,"|---:|---:|---:|---:|---:|---:|")
for r in peaks
    lm=nearest_lsts(v(r,"target_alpha"))
    println(io,"|$(v(r,"alpha_f"))|$(v(r,"kappa"))|$(v(r,"Gf"))|$(v(r,"target_alpha"))|$(lm.alpha)|$(lm.distance)|")
end
println(io)
println(io,"这些 alpha_f 粗扫描点的最大值如果位于 κ=0.02 或 κ=0.60 边界，只能说明扫描边界受限，不能拟合成内部匹配曲线。nearest spatial LST alpha 仅用于谱距离诊断，未将其标记为 Type I/II，也未将距离直接解释为共振。")
println(io)
println(io,"## 波数匹配判断")
println(io)
println(io,"当前数据需要同时比较 source_norm、G_vol 和 H_eff=sqrt(E_vol)/||S_vol||。如果 G_vol 与 source_norm 同步变化，增强主要来自粗糙度源项；如果 H_eff 在某个内部区域增强，才说明目标 LNS 传播/响应对强迫的相对放大发生变化。")
println(io,"改变 alpha_f 后必须使用每个 alpha_f 自己的低涡量入射根；本次远场根分析显示该根的 Im(kz) 从 alpha_f=0.20 的约 -0.0493 变化到 alpha_f=0.45 的约 -0.1934，不能固定使用 -0.1491i。")
println(io)
println(io,"## 收敛性检查")
println(io)
println(io,"|alpha_f|kappa|branch|G_full min/max over N,zmax|relative range|max qf derivative mismatch|")
println(io,"|---:|---:|---|---:|---:|---:|")
for (a,k,b) in [(0.38482,0.005,"plus"),(0.38482,0.4,"plus"),(0.38482,0.2,"minus"),(0.2,0.02,"plus"),(0.45,0.02,"plus")]
    rr=[r for r in conv if v(r,"alpha_f")==a && v(r,"kappa")==k && r["branch"]==b]
    gs=[v(r,"Gf") for r in rr]; spread=(maximum(gs)-minimum(gs))/maximum(gs); qd=maximum(v(r,"qf_deriv") for r in rr)
    println(io,"|$a|$k|$b|[$(minimum(gs)), $(maximum(gs))]|$spread|$qd|")
end
println(io)
println(io,"zmax 变化产生的约 (10^{-3}) 级 G_full 变化大于同一 zmax 下 N=89→129 的差异，说明当前代表性点的主要未确定性仍来自有限二阶外缘截断，而不是 Chebyshev 阶数。低 κ 壁面梯度在大 zmax 下仍需谨慎解释。")
println(io)
println(io,"## 固定目标波数实验")
println(io)
println(io,"固定负分支目标 alpha_hf=0.18482，令 kappa=alpha_f-0.18482：")
println(io,"|alpha_f|kappa|G_vol|G_full|H_eff|source norm|rho|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|")
for r in fixed
    println(io,"|$(v(r,"alpha_f"))|$(v(r,"kappa_abs"))|$(v(r,"Gv"))|$(v(r,"Gf"))|$(v(r,"H_eff"))|$(v(r,"source_norm"))|$(v(r,"rho"))|")
end
println(io,"固定目标实验使用 signed kappa=target-alpha_f，并显式断言 alpha_f+kappa_signed=target。alpha_f=0.38482、kappa_abs=0.20 的固定目标结果 G_full=$(v(first(r for r in fixed if isapprox(v(r,"alpha_f"),0.38482;atol=1e-12)),"Gf"))，与独立负分支扫描一致。固定目标下响应仍随 alpha_f 变化，说明负分支峰值不能只归因于目标 LNS 算子的固定波数响应；入射响应、粗糙度畸变和二阶源项本身也在系统变化。")
println(io)
println(io,"## 数值限制")
println(io)
println(io,"本阶段 alpha_f–kappa 扫描使用 N=89、zmax=30 的 Modal 粗扫描；低 κ 与正负分支使用 N=109、zmax=30。高响应和分支差异点还需要用 N=89,109,129 及 zmax=20,30,40 做收敛检查。二阶外缘仍是齐次四变量 Dirichlet，因此所有匹配结论都属于当前有限域模型。")
write(joinpath(OUTDIR,"wavenumber_matching_report.md"),String(take!(io)))
println("wrote report")

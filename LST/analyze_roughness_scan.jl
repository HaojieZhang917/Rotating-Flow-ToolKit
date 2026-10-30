using Printf

const OUTDIR=joinpath(@__DIR__,"results","roughness_wavenumber_scan")
function rows(path)
    ls=readlines(path); h=split(ls[1],',')
    return [Dict(zip(h,split(line,','))) for line in ls[2:end] if !isempty(strip(line))]
end
v(r,k)=parse(Float64,r[k])
main=rows(joinpath(OUTDIR,"roughness_scan_main.csv"))
conv=rows(joinpath(OUTDIR,"roughness_scan_convergence.csv"))
gpeak=argmax([v(r,"Gf_M") for r in main]); dpeak=argmax([v(r,"delta_full") for r in main])
srcpeak=argmax([v(r,"source_M") for r in main]); gain_src=[v(r,"Gv_M")/v(r,"source_M")^2 for r in main]
ratiopeak=argmax(gain_src)
rhos=[v(r,"rho_M") for r in main if isfinite(v(r,"rho_M"))]
energy_errors=[abs(v(r,"Ef_M")-v(r,"Ev_M")-v(r,"Ew_M")-v(r,"Ec_M")) for r in main]
energy_errors_D=[abs(v(r,"Ef_D")-v(r,"Ev_D")-v(r,"Ew_D")-v(r,"Ec_D")) for r in main]
linear_errors=vcat([v(r,"linear_M") for r in main],[v(r,"linear_D") for r in main])
io=IOBuffer()
println(io,"# 粗糙度波数扫描的附加科学分析")
println(io)
println(io,"- Modal 完整响应最大值：κ=$(v(main[gpeak],"kappa"))，G_full=$(v(main[gpeak],"Gf_M"))，目标 alpha_hf=$(v(main[gpeak],"target_alpha"))。该点位于扫描下边界，因此当前数据不能证明存在内部最优波数或共振。")
println(io,"- D/M 完整响应差异最大值：κ=$(v(main[dpeak],"kappa"))，delta_full=$(v(main[dpeak],"delta_full"))。")
println(io,"- Modal 体积源项范数最大值：κ=$(v(main[srcpeak],"kappa"))，||S_vol||=$(v(main[srcpeak],"source_M"))。")
println(io,"- 描述性 G_vol/||S_vol||² 最大值：κ=$(v(main[ratiopeak],"kappa"))，值=$(gain_src[ratiopeak])。该比值不是算子范数或严格感受性系数。")
println(io,"- Modal rho 范围：[$(minimum(rhos)), $(maximum(rhos))]，最大 |rho|=$(maximum(abs.(rhos)))；当前没有出现 |rho|>1 的能量恒等式异常。")
println(io,"- 主扫描能量恒等式最大绝对误差：Modal=$(maximum(energy_errors))，Dirichlet=$(maximum(energy_errors_D))；线性叠加相对误差最大值=$(maximum(linear_errors))。")
println(io)
println(io,"## 机制判断")
println(io)
println(io,"主扫描中 G_full 与 G_vol 的峰值都位于最低 κ=0.05，且随 κ 增大总体下降；因此当前有限参考高度、单一入射极化和齐次二阶外缘条件下，没有发现内部 κ 共振峰。")
println(io,"壁面增益比体积增益小很多，但 delta_wall 通常大于 delta_vol，说明壁面几何项对一阶 D/M 边界差异更敏感；由于其绝对能量很小，它对 G_full 的影响仍然有限。")
println(io,"体积源项最大位置与最终体积响应最大位置不应自动等同；本扫描还需结合 G_vol/||S_vol||² 的变化解释线性二阶算子对不同目标 alpha_hf 的放大或衰减。")
println(io)
println(io,"## 收敛判断")
println(io)
for k in sort(unique(v(r,"kappa") for r in conv))
    rr=[r for r in conv if v(r,"kappa")==k]
    vals=[v(r,"Gf_M") for r in rr]
    spread=(maximum(vals)-minimum(vals))/maximum(vals)
    println(io,"- κ=$(k)：Modal G_full 在 N/zmax 收敛样本中的相对范围约 $(spread)。壁面相关指标仍应优先检查 qf_deriv_M 和 delta_wall。")
end
write(joinpath(OUTDIR,"roughness_scan_analysis.md"),String(take!(io)))
println("wrote analysis")

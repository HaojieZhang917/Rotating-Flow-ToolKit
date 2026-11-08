#=
validate_turbulent_baseflow.jl
==============================

Checks of `cavity_baseflow_turb` (CavityLST.jl): Batchelor similarity equations with a stator-side
mixing-length eddy viscosity (Michel tanh mixing length + van Driest damping + Klebanoff intermittency),
laminar rotor layer and core, local similarity at radius R (delta units, local Omega r^2/nu = R^2).

 T1  nu_t = 0 reproduces the laminar solution;
 T2  grid convergence N = 200, 300, 400 (Re_h = 1000, R = 330);
 T3  R dependence at Re_h = 1000 (R = 100 ... 600);
 T4  sensitivity to the model constants (lam, A+) at Re_h = 1000, R = 330;
 T5  ORDER-OF-MAGNITUDE reference only (enclosed cavity with shroud and hub, developing non-similar layers;
     not a calibration or validation of the local-similarity model): DNS of Poncet & Randriamampianina (2006; C. R. Mecanique 2005), G = 18.32,
     a/b = 0.152, Re = 9.5e4 -> b/h = 21.6, Re_h = 203.6; r* = 0.56, 0.68, 0.80 -> R = 193, 225, 256.
     DNS: K = 0.375-0.418 (0.44 <= r* <= 0.8), v_tau/(Omega r) = 0.0341, 0.0335, 0.0355,
     delta_B+ ~ 47 at r* = 0.68 (z+ = 17 <-> z/delta_B = 0.365);
 T6  laminar K(Re_h) near Re_h = 204 (is the laminar reference smooth there?).
Output: results/turbulent_baseflow/{report.md, profiles.csv}
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
using LinearAlgebra, Printf
BLAS.set_num_threads(min(8, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "turbulent_baseflow"); mkpath(OUT)
io = open(joinpath(OUT, "report.md"), "w")
say(s) = (print(s); print(io, s); flush(io))
row(r) = @sprintf("%.4f | %.4f | %.4f | %.4f | %.4f | %.2f | %.1f | %.2f | %.4f", r.K, r.utau, r.delta_s, r.delta_s_peak, r.delta_r,
                  r.delta_s_peak/r.delta_r, r.delta_s_plus, maximum(r.nut), r.base.Pi)
hdr = "K | u_τ/(Ωr) | δ_s(模型外尺度) | δ_s(峰值) | δ_r | δ_s(峰值)/δ_r | δ_s⁺ | max ν_t/ν | Π"

say("# 湍流静盘层基流（局部相似 + 混合长涡黏性）验证\n\n")
say("模型：ν_t/ν = (ℓ/h)² R√Re_h √(F′²+G′²) γ，ℓ/h = λδ_s tanh(κyD/(λδ_s))，D = 1−exp(−y⁺/A⁺)，γ = 1/(1+5.5(y/δ_s)⁶)，y = 1−z；默认 λ = 0.085、A⁺ = 26、κ = 0.41（平衡湍流边界层的标准常数，未针对本问题调参）；δ_s 取静盘侧 G 超调峰值到静盘的距离（edge = :peak，理由：层流/弱湍流 Bödewadt 层带超调，0.99K 首次到达点只在层的内半部；:first99 作为敏感性对照）。δ_r：G 降到 K+0.01(1−K) 处到转盘的距离（间隙单位）。\n\n")

# T1
bl = cavity_baseflow(1000.0; N=300)
b0 = solve_cavity_baseflow(1000.0; N=300, guess=bl, nut=zeros(301))
say(@sprintf("## T1 ν_t = 0\n\n|ΔG|∞ = %.1e，|ΔΠ| = %.1e（相对层流解）\n\n", maximum(abs.(b0.G - bl.G)), abs(b0.Pi - bl.Pi)))

# T2
say("## T2 网格收敛（Re_h = 1000，R = 330）\n\n| N | $hdr |\n|---|---|---|---|---|---|---|---|---|---|\n")
for N in (200, 300, 400)
    r = cavity_baseflow_turb(1000.0, 330.0; N=N)
    say("| $N | " * replace(row(r), " | " => " | ") * " |\n")
end
say("\n")

# T3
say("## T3 R 的影响（Re_h = 1000；层流 K = 0.3147）\n\n| R | Ωr²/ν | $hdr |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
prof = Dict{Float64,Any}()
for R in (100.0, 200.0, 330.0, 450.0, 600.0)
    r = cavity_baseflow_turb(1000.0, R); prof[R] = r
    say(@sprintf("| %.0f | %.2e | ", R, R^2) * row(r) * " |\n")
end
say("\n")

# T4
say("## T4 模型常数敏感性（Re_h = 1000，R = 330）\n\n| λ | A⁺ | δ_s 定义 | $hdr |\n|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (lam, Ap, ed) in ((0.07, 26.0, :first99), (0.085, 26.0, :first99), (0.10, 26.0, :first99), (0.085, 20.0, :first99),
                      (0.085, 32.0, :first99), (0.085, 26.0, :peak))
    r = cavity_baseflow_turb(1000.0, 330.0; lam=lam, Aplus=Ap, edge=ed)
    say(@sprintf("| %.3f | %.0f | %s | ", lam, Ap, ed) * row(r) * " |\n")
end
say("\n")

# T5
say("## T5 量级参考：Poncet & Randriamampianina DNS（Re_h = 203.6）\n\n注意：该 DNS 是带外壁和内轴的封闭腔，静盘层从外壁向内发展、不满足局部相似，K 受约束效应和径向发展历史影响。本模型是无限大双盘的局部相似解，二者不可直接比较；此表只作壁面量（u_τ、δ⁺）的量级参考，不用于标定，K 的对比没有意义。\n\n")
lam204 = cavity_baseflow(203.6; N=300); e = stator_eddy_viscosity(lam204, 200.0; lam=0.0, edge=:first99)
say(@sprintf("层流相似解：K = %.4f，δ_s(0.99K) = %.4f\n\n| r* | R | δ_s 定义 | 模型 K | DNS K | 模型 u_τ/(Ωr) | DNS u_τ/(Ωr) | 模型 δ_s⁺ | δ_s(峰值)/δ_r |\n|---|---|---|---|---|---|---|---|---|\n", e.K, e.delta_s))
for (rs, R, ut) in ((0.56, 193.0, 0.0341), (0.68, 225.0, 0.0335), (0.80, 256.0, 0.0355))
    for ed in (:first99, :peak)
        r = cavity_baseflow_turb(203.6, R; edge=ed)
        say(@sprintf("| %.2f | %.0f | %s | %.4f | 0.375–0.418 | %.4f | %.4f | %.1f | %.2f |\n", rs, R, ed, r.K, r.utau, ut, r.delta_s_plus, r.delta_s_peak/r.delta_r))
    end
end
say("\nDNS δ_B⁺ ≈ 47（r* = 0.68，由 z⁺ = 17 ↔ z/δ_B = 0.365 换算）。\n\n")

# T6
say("## T6 层流 K(Re_h)\n\n| Re_h | K |\n|---|---|\n")
for Re in (150.0, 175.0, 203.6, 250.0, 300.0, 400.0, 600.0, 1000.0)
    b = cavity_baseflow(Re; N=300); say(@sprintf("| %.1f | %.4f |\n", Re, bary_interp(b.z, b.G, 0.5)))
end
close(io)

open(joinpath(OUT, "profiles.csv"), "w") do f
    println(f, "case,z,F,G,H,nut")
    for j in eachindex(bl.z); println(f, "laminar,$(bl.z[j]),$(bl.F[j]),$(bl.G[j]),$(bl.H[j]),0"); end
    for R in sort(collect(keys(prof))), j in eachindex(prof[R].base.z)
        r = prof[R]
        println(f, "R$(Int(R)),$(r.base.z[j]),$(r.base.F[j]),$(r.base.G[j]),$(r.base.H[j]),$(r.nut[j])")
    end
end
println("\nwrote ", OUT)

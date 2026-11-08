#=
validate_cavity_port.jl
=======================

Validation of CavityLST.jl against the manuscript (Zhang & Liu, rotor–stator, surface
mass flux) and the stored base flow TwoDisk/baseflow_Res1000.npz.

  B1  base flow at Re_h = 1000, a_s = 0: core G_c, Pi; pointwise comparison with the npz
      profiles (u0 = F, v0 = G, w0 = H, du0 = F', dv0 = G'); base-grid convergence.
  E1  grid-convergence case of Table I: Re_h = 2000, R = 330, beta = 0.0503, omega_bar = 8:
      alpha = 0.457945238 - 0.013718727i (N_z = 300).
  E2  neutral Type-II mode, Re_h = 1000, R = 100, n = -7, omega_bar = 8: alpha ≈ 0.272.
  E3  mode identification, R = 160, n = 16, omega = 0.05: 0.174 + 0.608i (rotor), 0.342 - 0.175i (stator).
  E4  stator convective critical point: R = 47.0, beta = -0.009, omega = 0.020, alpha ≈ 0.281.
  E5  rotor stationary lower neutral point for n = 30: R_l = 287.822.
  C1  local receptivity coefficient, omega_bar = 0, n = 30, c = 1, hr = 1/Re_h:
      |C_r(R_l = 287.822)| = 0.1436, maximum 0.2444 at R_p = 468.910.
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
include(joinpath(@__DIR__, "local_mode_tracking.jl"))
using LinearAlgebra, Printf
BLAS.set_num_threads(min(16, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "cavity_port_validation"); mkpath(OUT)
io = IOBuffer(); say(s) = (println(s); println(io, s); flush(stdout))
say("# CavityLST 移植验证（对照论文与 TwoDisk 原代码）\n")

# ---------------- B1 base flow
say("## B1 基流（Re_h=1000，a_s=0）\n")
b = cavity_baseflow(1000.0; N=300)
Gc = bary_interp(b.z, b.G, 0.5)
say(@sprintf("- Newton 残差 %.1e（%d 次迭代）；G_c=G(0.5)=%.5f（论文 0.313）；Π=%.5f（论文 ≈ -0.098）；-G_c²=%.5f。", b.residual, b.iterations, Gc, b.Pi, -Gc^2))
npz = let
    py = """
import numpy as np, sys
d=np.load(sys.argv[1]);
for k in ['x','u0','v0','w0','du0','dv0']: np.savetxt(sys.argv[2]+'_'+k+'.txt', d[k])
"""
    tmp = joinpath(OUT, "npz"); write(tmp*"_conv.py", py)
    run(`python3 -I $(tmp*"_conv.py") $(joinpath(@__DIR__, "..", "TwoDisk", "baseflow_Res1000.npz")) $tmp`)
    Dict(k => parse.(Float64, readlines(tmp*"_"*k*".txt")) for k in ["x", "u0", "v0", "w0", "du0", "dv0"])
end
let x = npz["x"], idx = 1:50:length(npz["x"])
    DF = b.D*b.F; DG = b.D*b.G
    e = Dict(k => maximum(abs(bary_interp(b.z, f, x[i]) - npz[nk][i]) for i in idx) for (k, f, nk) in
             (("F", b.F, "u0"), ("G", b.G, "v0"), ("H", b.H, "w0"), ("F'", DF, "du0"), ("G'", DG, "dv0")))
    say(@sprintf("- 与 baseflow_Res1000.npz 逐点最大差：F %.1e，G %.1e，H %.1e，F' %.1e，G' %.1e；F'(0)=%.5f（npz %.5f），G'(0)=%.5f（npz %.5f）。",
            e["F"], e["G"], e["H"], e["F'"], e["G'"], DF[1], npz["du0"][1], DG[1], npz["dv0"][1]))
end
for Nb in (200, 400)
    bb = cavity_baseflow(1000.0; N=Nb)
    say(@sprintf("- 基流网格 N=%d：|ΔG_c|=%.1e，|ΔΠ|=%.1e，|ΔF'(0)|=%.1e。", Nb, abs(bary_interp(bb.z, bb.G, 0.5) - Gc), abs(bb.Pi - b.Pi), abs((bb.D*bb.F)[1] - (b.D*b.F)[1])))
end
say("")

nearest_mode(op, target) = (sp = cavity_spectrum(op; shift=target); m = refine_mode(op.L0, op.L1, op.L2, sp[1].alpha); (m, sp[1].alpha))

# ---------------- E1
say("## 特征值对照\n")
say("|算例|参数|本文 α|论文 α|差|\n|---|---|---|---|---:|")
b2000 = cavity_baseflow(2000.0; N=400)
for N in (199, 299)
    op = cavity_pencil(b2000, 330.0, 0.0503, 8/330; N=N)
    m, _ = nearest_mode(op, 0.4579 - 0.0137im)
    ref = 0.457945238 - 0.013718727im
    say(@sprintf("|E1 网格收敛（N_z=%d）|Re_h=2000，R=330，β=0.0503，ω̄=8|%.9f%+.9fi|%.9f%+.9fi|%.1e|", N+1, real(m.alpha), imag(m.alpha), real(ref), imag(ref), abs(m.alpha - ref)))
end
let op = cavity_pencil(b, 100.0, -7/100, 8/100; N=199), (m, _) = nearest_mode(op, 0.272 + 0im)
    say(@sprintf("|E2 Type-II 中性|R=100，n=-7，ω̄=8|%.6f%+.6fi|0.272（中性）|%.1e|", real(m.alpha), imag(m.alpha), abs(m.alpha - 0.272)))
end
let op = cavity_pencil(b, 160.0, 16/160, 0.05; N=199)
    for ref in (0.174 + 0.608im, 0.342 - 0.175im)
        m, _ = nearest_mode(op, ref)
        f = full_state(op, m.vector); n = op.n; vel = abs.(f[1:n]) .+ abs.(f[n+1:2n]) .+ abs.(f[2n+1:3n])
        zpk = op.z[argmax(vel)]
        say(@sprintf("|E3 模态辨识|R=160，n=16，ω=0.05|%.4f%+.4fi（峰值 z=%.3f）|%.3f%+.3fi|%.1e|", real(m.alpha), imag(m.alpha), zpk, real(ref), imag(ref), abs(m.alpha - ref)))
    end
end
let op = cavity_pencil(b, 47.0, -0.009, 0.020; N=199), (m, _) = nearest_mode(op, 0.281 + 0im)
    say(@sprintf("|E4 静盘临界点|R=47.0，β=-0.009，ω=0.020|%.5f%+.5fi|0.281（中性，表中参数取三位有效数字）|%.1e|", real(m.alpha), imag(m.alpha), abs(m.alpha - 0.281)))
end
let op = cavity_pencil(b, 287.822, 30/287.822, 0.0; N=199), (m, _) = nearest_mode(op, 0.40 + 0im)
    say(@sprintf("|E5 n=30 下中性点|R=287.822，ω̄=0|%.6f%+.2ei|Im α=0|%.1e|", real(m.alpha), imag(m.alpha), abs(imag(m.alpha))))
end
say("")

# ---------------- C1 receptivity
say("## C1 局部感受性系数（ω̄=0，n=30，c=1，h_r=1/Re_h）\n")
function Cr(b, R, guess; N=199)
    op = cavity_pencil(b, R, 30/R, 0.0; N=N)
    m = refine_mode(op.L0, op.L1, op.L2, guess)
    a = refine_mode(transpose(op.L0), transpose(op.L1), transpose(op.L2), m.alpha)
    f = full_state(op, m.vector); n = op.n
    x = m.vector ./ f[n+argmax(abs.(f[n+1:2n]))]                 # max|v| = 1, as in the original code
    Λ = residue_coefficient(op, m.alpha, x, a.vector, wall_roughness_rhs(op, m.alpha))
    return abs(Λ), m.alpha
end
op0 = cavity_pencil(b, 287.822, 30/287.822, 0.0; N=199)
a0 = nearest_mode(op0, 0.40 + 0im)[1].alpha
c_l, al = Cr(b, 287.822, a0)
say(@sprintf("- R_l=287.822：|C_r|=%.4f（论文 0.1436），α=%.6f%+.2ei。", c_l, real(al), imag(al)))
Rs = collect(287.822:4.0:520.0); vals = Float64[]; a = a0
for R in Rs
    global a
    c, a = Cr(b, R, a); push!(vals, c)
end
i = argmax(vals)
# refine the maximum with a parabola through the three neighbours
Rp = if 1 < i < length(Rs)
    x1, x2, x3 = Rs[i-1:i+1]; y1, y2, y3 = vals[i-1:i+1]
    x2 - 0.5*((x2-x1)^2*(y2-y3) - (x2-x3)^2*(y2-y1))/((x2-x1)*(y2-y3) - (x2-x3)*(y2-y1))
else
    Rs[i]
end
cp, _ = Cr(b, Rp, a)
say(@sprintf("- 最大值：R_p≈%.2f，|C_r|=%.4f（论文 R_p=468.910，0.2444）。", Rp, cp))
for N in (149, 249)
    c, _ = Cr(b, 287.822, a0; N=N)
    say(@sprintf("- N=%d 时 |C_r(R_l)|=%.6f（与 N=199 相差 %.1e）。", N, c, abs(c - c_l)/c_l))
end
open(joinpath(OUT, "Cr_n30_omega0.csv"), "w") do f
    println(f, "R,abs_Cr"); for (R, v) in zip(Rs, vals); println(f, "$R,$v"); end
end
write(joinpath(OUT, "cavity_port_validation.md"), String(take!(io)))
println("wrote ", OUT)

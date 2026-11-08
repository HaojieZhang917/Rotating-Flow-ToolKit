#=
analyze_gap_filter.jl
=====================

Verification of the "gap filter" theory (results/gap_filter/gap_filter_theory.md) with the full
rotor–stator operator of CavityLST.jl (rotor frame, local parallel, P_N-P_{N-2}).

Theory (inviscid rotating core, frame co-rotating with the core, local Cartesian):
  P'' = m^2 P,  m^2 = k_h^2 sigma^2/(sigma^2 - f^2),  k_h^2 = alpha^2 + beta^2,
  sigma = omega + beta (1 - K)   (Doppler shift rotor frame -> core frame),  f = 2K/R  (local units),
  rigid lid (w = 0) at the rotor-layer edge  ->  T = u_h(rotor edge)/u_h(stator edge) = 1/cosh(m Dzeta),
  Dzeta = (z_s - z_r) sqrt(Re_h)  (delta units).

 A  exact roots of the constant-coefficient core operator (mid-core F, G, H; viscosity, axial flow,
    curvature included) versus the inviscid +-m;
 B  vertical exponents of the full numerical solution in the core (two-exponential linear-prediction
    fit of u on 0.25 <= z <= 0.75) versus A and the inviscid formula;
 C  end-to-end transfer |u_h(z_r)|/|u_h(z_s)| versus 1/|cosh(m Dzeta)|;
 D  Re_h scan (1000, 2000, 4000) of the transfer at fixed (alpha, n, c);
 E  direct-path coefficient Lambda_dir of the Type-I target (R = 330, n = 16, omega_bar = 8) versus Re_h,
    compared with exp(-m Dzeta);
 F  identification of the near-resonance peaks of the incident field (sigma/f, |m| Dzeta, and the
    location of the nearest eigenmode of the incident pencil).
Forcing: stator-layer Gaussian e_r at z_s0 = 1 - 0.045 (sigma_z = 0.02), omega = -c beta.
Output: results/gap_filter/{report.md, partB.csv, partC.csv, partD.csv, partE.csv}
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
include(joinpath(@__DIR__, "local_mode_tracking.jl"))
using LinearAlgebra, Printf
BLAS.set_num_threads(min(8, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "gap_filter"); mkpath(OUT)
io = open(joinpath(OUT, "report.md"), "w")
say(s) = (print(s); print(io, s); flush(io))

# ------------------------------------------------------------------ helpers
pmat(op, a) = op.raw[1] + a*op.raw[2] + a^2*op.raw[3]
function forced(b, R, alpha, n, c; N=199)
    beta = n/R; op = cavity_pencil(b, R, beta, -c*beta; N=N, pressure=:PN2)
    f = zeros(ComplexF64, 4op.n); f[1:op.n] .= exp.(-((op.z .- 0.955).^2) ./ (2*0.02^2))
    q = zeros(ComplexF64, 4op.n); P = pmat(op, alpha); q[op.keep] .= P[op.keep, op.keep] \ f[op.keep]
    return op, q
end
comp(op, q, j) = q[(j-1)*op.n+1:j*op.n]
uh(op, q, z) = sqrt(abs2(bary_interp(op.z, comp(op, q, 1), z)) + abs2(bary_interp(op.z, comp(op, q, 2), z)))

"""Layer edges: rotor = first G minimum above the rotor, stator = G maximum below the stator (overshoot)."""
function edges(b)
    zq = range(0, 0.5, length=5001); Gq = [bary_interp(b.z, b.G, t) for t in zq]
    k = findfirst(i -> Gq[i] <= Gq[i-1] && Gq[i] <= Gq[i+1], 2:length(Gq)-1)
    yq = range(0, 0.5, length=5001); ys = yq[argmax([bary_interp(b.z, b.G, 1 - t) for t in yq])]
    return zq[k+1], 1 - ys
end

"""Inviscid vertical exponent m (principal branch, Re m >= 0) in delta units."""
function m_inviscid(alpha, beta, omega, K, R)
    sig = omega + beta*(1 - K); f = 2K/R; kh2 = alpha^2 + beta^2
    return sqrt(complex(kh2*sig^2/(sig^2 - f^2))), sig, f
end

"""
Exact roots lambda (delta units, q ~ exp(lambda zeta)) of the constant-coefficient core operator at the
mid-core values (F, G-1, H) at z = 1/2 (rotor frame, same terms as cavity_pencil with D -> sqrt(Re_h) lambda).
"""
function core_roots(b, R, alpha, beta, omega)
    Fc = bary_interp(b.z, b.F, 0.5); Gc = bary_interp(b.z, b.G, 0.5) - 1; Hc = bary_interp(b.z, b.H, 0.5)
    d0 = -im*omega + im*beta*Gc + im*alpha*Fc + (alpha^2 + beta^2)/R + Fc/R
    # M(l) = M0 + l M1 + l^2 M2 acting on (u, v, w, p)
    M0 = ComplexF64[d0 -2(Gc+1)/R 0 im*alpha; 2(Gc+1)/R d0 0 im*beta; 0 0 d0 0; 1/R+im*alpha im*beta 0 0]
    M1 = ComplexF64[Hc/R 0 0 0; 0 Hc/R 0 0; 0 0 Hc/R 1; 0 0 1 0]
    M2 = ComplexF64[-1/R 0 0 0; 0 -1/R 0 0; 0 0 -1/R 0; 0 0 0 0]
    Z = zeros(ComplexF64, 4, 4); Id = Matrix{ComplexF64}(I, 4, 4)
    E = eigen([-M1 -M0; Id Z], [M2 Z; Z Id]).values
    return sort([l for l in E if isfinite(l) && abs(l) < 1e3], by=abs)
end

"""Two-exponential linear-prediction (Prony) fit of a field on uniform zeta samples; returns the two exponents."""
function prony2(zeta, u)
    h = zeta[2] - zeta[1]; M = length(u)
    A = hcat(u[2:M-1], u[1:M-2]); y = u[3:M]
    c = A \ y; r = roots2(c[1], c[2])
    lam = log.(complex.(r)) ./ h
    res = norm(A*c - y)/norm(y)
    return lam, res
end
roots2(c1, c0) = (d = sqrt(complex(c1^2 + 4c0)); [(c1 + d)/2, (c1 - d)/2])

Re0, R0 = 1000.0, 330.0
b0 = cavity_baseflow(Re0; N=300); K0 = bary_interp(b0.z, b0.G, 0.5); zr0, zs0 = edges(b0)
say("# 间隙滤波理论的数值验证\n\nRe_h = $(Re0)，R = $(R0)，K = $(round(K0, digits=4))；转盘层外缘 z_r = $(round(zr0, digits=4))（G 的第一个极小），静盘层外缘 z_s = $(round(zs0, digits=4))（G 的超调极大）。f = 2K/R = $(round(2K0/R0, sigdigits=4))。\n\n")

# ------------------------------------------------------------------ A
say("## A 核心常系数问题的精确根 vs 无黏公式\n\n| α | n | c | σ/f | 无黏 m | 精确根（|λ| 最小的两个） |\n|---|---|---|---|---|---|\n")
for (a, n, c) in ((0.45, 16, 8/16*0+0.0), (0.05, -10, 0.8), (0.05, -10, 0.7), (0.1, -20, 0.69), (0.1, 0, 0.8), (0.3, 10, 0.85))
    beta = n/R0; om = n == 16 ? 8/R0 : -c*beta
    m, sig, f = m_inviscid(a, beta, om, K0, R0); lr = core_roots(b0, R0, a, beta, om)
    say(@sprintf("| %.2f | %d | %s | %.3g | %.4f%+.4fi | %s |\n", a, n, n == 16 ? "ω̄=8" : @sprintf("%.2f", c), real(sig/f), real(m), imag(m),
                 join([@sprintf("%.4f%+.4fi", real(l), imag(l)) for l in lr[1:2]], "，")))
end
say("\n（精确根包含黏性、核心轴向流 H 和曲率；另外几个根是 |λ| ≫ 1 的黏性根，只存在于薄层内。）\n\n")

# ------------------------------------------------------------------ B and C
say("## B/C 全数值解：核心竖向指数与端到端传递比\n\n")
rowsB = String[]
AF = vcat(0.01, 0.02:0.02:0.6); NFs = (-20, -10, 0, 10, 20); Cs = (0.7, 0.75, 0.8, 0.9, 1.0)
zfit = range(0.25, 0.75, length=61); zeta = collect(zfit) .* sqrt(Re0)
for n in NFs, c in Cs, a in AF
    n == 0 && c != 0.8 && continue                       # sigma = 0 for n = 0, c irrelevant
    beta = n/R0; om = -c*beta
    op, q = forced(b0, R0, a, n, c)
    u = [bary_interp(op.z, comp(op, q, 1), t) for t in zfit]
    v = [bary_interp(op.z, comp(op, q, 2), t) for t in zfit]
    w = norm(u) > norm(v) ? u : v
    lam, res = prony2(zeta, w)
    m, sig, f = m_inviscid(a, beta, om, K0, R0); lr = core_roots(b0, R0, a, beta, om)[1:2]
    Tn = uh(op, q, zr0)/uh(op, q, zs0); Dz = (zs0 - zr0)*sqrt(Re0)
    Tt = 1/abs(cosh(m*Dz))
    push!(rowsB, join((a, n, c, real(sig/f), real(m), imag(m), real(lr[1]), imag(lr[1]), real(lr[2]), imag(lr[2]),
                       real(lam[1]), imag(lam[1]), real(lam[2]), imag(lam[2]), res, Tn, Tt, Dz), ","))
end
open(joinpath(OUT, "partB.csv"), "w") do f
    println(f, "alpha,n,c,sig_over_f,m_r,m_i,ex1_r,ex1_i,ex2_r,ex2_i,fit1_r,fit1_i,fit2_r,fit2_i,fit_res,T_num,T_theory,Dzeta")
    foreach(r -> println(f, r), rowsB)
end
# summary of B: relative error of the fitted decay rate vs |Re m| for clearly evanescent cases
let errs = Float64[], erre = Float64[], cnt = 0
    for r in rowsB
        s = parse.(Float64, split(r, ",")); sf = s[4]; mr = s[5]; fitr = maximum(abs.((s[11], s[13])))
        exr = maximum(abs.((s[7], s[9]))); res = s[15]
        (abs(sf) > 1.2 && mr > 0.02 && res < 1e-3) || continue
        cnt += 1; push!(errs, abs(fitr - mr)/mr); push!(erre, abs(fitr - exr)/exr)
    end
    say(@sprintf("衰减区（|σ/f| > 1.2、m > 0.02、拟合残差 < 1e-3）共 %d 个点：拟合指数 vs 无黏 m 的相对误差 中位数 %.2e、最大 %.2e；vs 精确根 中位数 %.2e、最大 %.2e。\n\n",
                 cnt, sort(errs)[max(1, cnt ÷ 2)], maximum(errs), sort(erre)[max(1, cnt ÷ 2)], maximum(erre)))
end

# ------------------------------------------------------------------ D
say("## D Re_h 扫描：端到端传递比\n\n| Re_h | K | Δζ | α | n | c | m | ln T_num | ln T_theory |\n|---|---|---|---|---|---|---|---|---|\n")
rowsD = String[]
for Re in (1000.0, 2000.0, 4000.0)
    b = cavity_baseflow(Re; N=400); K = bary_interp(b.z, b.G, 0.5); zr, zs = edges(b)
    Dz = (zs - zr)*sqrt(Re)
    for (a, n, c) in ((0.1, 10, 0.85), (0.2, -10, 1.0), (0.3, 0, 0.8), (0.3, 20, 0.9), (0.45, 16, -0.5))
        beta = n/R0; om = -c*beta
        op, q = forced(b, R0, a, n, c; N=249)
        m, = m_inviscid(a, beta, om, K, R0)
        Tn = uh(op, q, zr)/uh(op, q, zs); Tt = 1/abs(cosh(m*Dz))
        say(@sprintf("| %.0f | %.4f | %.1f | %.2f | %d | %.2f | %.4f | %.2f | %.2f |\n", Re, K, Dz, a, n, c, real(m), log(Tn), log(Tt)))
        push!(rowsD, join((Re, K, Dz, a, n, c, real(m), log(Tn), log(Tt)), ","))
    end
end
open(joinpath(OUT, "partD.csv"), "w") do f
    println(f, "Re_h,K,Dzeta,alpha,n,c,m,lnT_num,lnT_theory"); foreach(r -> println(f, r), rowsD)
end
say("\n")

# ------------------------------------------------------------------ E
say("## E 直接路径 Λ_dir 与 Re_h（Type-I：R = 330，n_t = 16，ω̄ = 8）\n\n| Re_h | α_n | Δζ | m | ln|Λ_dir| | −mΔζ | ln|Λ_dir| + mΔζ |\n|---|---|---|---|---|---|---|\n")
rowsE = String[]; let a0 = 0.448780 - 0.014605im
    for Re in (1000.0, 1500.0, 2000.0, 3000.0, 4000.0)
        b = cavity_baseflow(Re; N=400); K = bary_interp(b.z, b.G, 0.5); zr, zs = edges(b)
        op = cavity_pencil(b, R0, 16/R0, 8/R0; N=249, pressure=:PN2)
        md = refine_mode(op.L0, op.L1, op.L2, a0); an = md.alpha; a0 = an
        y = refine_mode(transpose(op.L0), transpose(op.L1), transpose(op.L2), an).vector
        xf = full_state(op, md.vector); n = op.n
        zz = range(0, 0.3, length=3001); vz(t) = bary_interp(op.z, xf[n+1:2n], t)
        zst = 0.0478*sqrt(1000/Re)                          # same height in delta units as at Re_h = 1000
        x = md.vector ./ vz(zst)
        f = zeros(ComplexF64, 4n); f[1:n] .= exp.(-((op.z .- (1 - 0.045*sqrt(1000/Re))).^2) ./ (2*(0.02*sqrt(1000/Re))^2))
        L = residue_coefficient(op, an, x, y, f[op.keep])
        Dz = (zs - zr)*sqrt(Re); m, = m_inviscid(real(an), 16/R0, 8/R0, K, R0)
        say(@sprintf("| %.0f | %.5f%+.5fi | %.1f | %.4f | %.2f | %.2f | %.2f |\n", Re, real(an), imag(an), Dz, real(m), log(abs(L)), -real(m)*Dz, log(abs(L)) + real(m)*Dz))
        push!(rowsE, join((Re, real(an), imag(an), Dz, real(m), log(abs(L))), ","))
    end
end
open(joinpath(OUT, "partE.csv"), "w") do f
    println(f, "Re_h,alpha_n_r,alpha_n_i,Dzeta,m,lnLdir"); foreach(r -> println(f, r), rowsE)
end
say("\n源项的位置和宽度按 δ 单位固定（距静盘 0.045h、宽 0.02h，换算自 Re_h = 1000），以免源项本身随 Re_h 变化。\n\n")

# ------------------------------------------------------------------ F
say("## F 入射场近共振峰的识别（Re_h = 1000，R = 330）\n\n")
say("| n | c | 峰值 α_f | 峰值 T_num | σ/f | m 或 |m| | 最近特征值 α_j | 该模态 |v| 的峰值位置 z |\n|---|---|---|---|---|---|---|---|\n")
for n in NFs, c in Cs
    n == 0 && c != 0.8 && continue
    sel = [parse.(Float64, split(r, ",")) for r in rowsB if (s = parse.(Float64, split(r, ",")); s[2] == n && s[3] == c)]
    Tn = [s[16] for s in sel]; as = [s[1] for s in sel]
    # local maxima of T_num exceeding the theory by a factor > 3
    for k in 2:length(Tn)-1
        (Tn[k] > Tn[k-1] && Tn[k] > Tn[k+1] && Tn[k] > 3*sel[k][17]) || continue
        beta = n/R0; om = -c*beta
        op = cavity_pencil(b0, R0, beta, om; N=199, pressure=:PN2)
        sp = cavity_spectrum(op; shift=as[k] + 0im, maxabs=3.0)
        j = argmin([abs(s.alpha - as[k]) for s in sp]); xj = full_state(op, sp[j].vector)
        vj = abs.(xj[op.n+1:2op.n]); zpk = op.z[argmax(vj)]
        m, sig, f = m_inviscid(as[k], beta, om, K0, R0)
        say(@sprintf("| %d | %.2f | %.2f | %.3g | %.2f | %.3f | %.4f%+.4fi | %.3f |\n", n, c, as[k], Tn[k], real(sig/f), abs(m), real(sp[j].alpha), imag(sp[j].alpha), zpk))
    end
end
close(io)
println("\nwrote ", OUT)

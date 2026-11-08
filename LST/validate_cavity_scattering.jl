#=
validate_cavity_scattering.jl
=============================

Validation of the second-order (stator-forcing x rotor-roughness) scattering problem in
CavityLST.jl (P_N-P_{N-2} pressure throughout), Re_h = 1000, a_s = 0, R = 330.

Target  : rotor Type-I mode, n_t = 16, omega_bar_t = 8 (rotor frame), spatially growing alpha_n.
Incident: stator-layer forcing f_s = e_r exp[-(z-0.955)^2/(2*0.02^2)] at (alpha_f, n_f = -10,
          omega_f = omega_t) -> convection speed c = -omega_bar_t/n_f = 0.8 (in [1-G_c, 1]).
Roughness: rotor-steady Fourier wave (kappa = alpha_n - alpha_f, n_h = n_t - n_f = 26), unit amplitude.
Chain   : q_f = P_f \ f_s;  q_h = P_h \ wall(-F'(0), -G'(0), 0);
          P_t(alpha_t) q_hf = -(B(q_h,q_f) + B(q_f,q_h)) + wall q_hf(0) = -dq_f/dz(0);
          Lambda_hf = y^T b(alpha_n)/Q (complex-alpha residue, kappa complex).
Competing direct path: Lambda_dir = y^T f_s(alpha_n, n_t, omega_t)/Q.

Checks
  V1  bilinear form with a = base flow (F, G-1, H/R) == convective part of the operator (rotor frame,
      Coriolis 2/R separated), both pressure formulations;
  V2  complex detuning of the full second-order problem: delta q_hf(alpha_n + delta) -> x Lambda_hf
      (q_h recomputed at kappa = alpha_n + delta - alpha_f), error must decrease ~ |delta|;
  V3  grid convergence of Lambda_hf (volume / wall / total) and Lambda_dir, N = 99 ... 249.
x is normalised by v(z*) = 1 at a fixed height z* (max |v| of the N = 199 mode).
Output: results/cavity_scattering_validation/{report.md, lambda_convergence.csv}
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
include(joinpath(@__DIR__, "local_mode_tracking.jl"))
using LinearAlgebra, Printf, Random
BLAS.set_num_threads(min(16, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "cavity_scattering_validation"); mkpath(OUT)
io = open(joinpath(OUT, "report.md"), "w")
say(s) = (print(s); print(io, s); flush(io))

const Re_h, R = 1000.0, 330.0
const n_t, wbar = 16, 8.0
const n_f = -10; const n_h = n_t - n_f
const ZS, SIG = 0.955, 0.02
const ALPHA0 = 0.448780 - 0.014605im
b = cavity_baseflow(Re_h; N=300)
bt, bh, bf = n_t/R, n_h/R, n_f/R; om = wbar/R

say("# 二阶散射（静盘强迫 × 转盘粗糙元）验证\n\n")
say(@sprintf("Re_h = %.0f，R = %.0f；目标 Type-I：n_t = %d，ω̄_t = %.0f；入射 n_f = %d（c = %.3f），粗糙元 n_h = %d，ω_h = 0。全部使用 P_N–P_{N−2} 压力格式。\n\n",
             Re_h, R, n_t, wbar, n_f, -wbar/n_f, n_h))

# ------------------------------------------------------------------ V1
say("## V1 双线性项 vs 算子对流部分\n\n")
for pr in (:PN, :PN2)
    α, β = 0.37, 0.05
    op = cavity_pencil(b, R, β, 0.0; N=80, pressure=pr); n = op.n; D = op.D
    Random.seed!(3); q = randn(ComplexF64, 4n); q[3n+1:4n] .= 0          # p = 0: no pressure terms
    Lq = (op.raw[1] + α*op.raw[2] + α^2*op.raw[3])*q
    visc(x) = (α^2 + β^2)/R .* x .- (1/(R*Re_h)).*(op.D*(op.D*x))
    u, v, w = q[1:n], q[n+1:2n], q[2n+1:3n]
    cu = Lq[1:n] .- visc(u) .+ (2/R).*v                                  # remove viscosity and Coriolis
    cv = Lq[n+1:2n] .- visc(v) .- (2/R).*u
    cw = Lq[2n+1:3n] .- visc(w)
    F, Gl, H = op.F, op.G, op.H
    U = vcat(F, Gl .- 1, H ./ R, zeros(n))
    N1 = cavity_bilinear(U, q, α, β, op); N2 = cavity_bilinear(q, U, 0.0, 0.0, op)
    bu = N1[1] .+ N2[1] .+ (F/R).*u                                       # + radial derivative of r F
    bv = N1[2] .+ N2[2] .+ ((Gl .- 1)/R).*u                               # + radial derivative of r (G-1)
    bw = N1[3] .+ N2[3]
    say(@sprintf("- pressure = %s：相对误差 %.2e\n", pr, norm([bu - cu; bv - cv; bw - cw])/norm([cu; cv; cw])))
end
say("\n")

# ------------------------------------------------------------------ building blocks
function setup(N)
    opt = cavity_pencil(b, R, bt, om; N=N, pressure=:PN2)
    m = refine_mode(opt.L0, opt.L1, opt.L2, ALPHA0)
    an = m.alpha
    y = refine_mode(transpose(opt.L0), transpose(opt.L1), transpose(opt.L2), an).vector
    return opt, an, m.vector, y
end
function stator_forcing(op)
    f = zeros(ComplexF64, 4op.n); f[1:op.n] .= exp.(-((op.z .- ZS).^2)./(2SIG^2)); f
end
solve(op, α, rhs) = (q = zeros(ComplexF64, 4op.n); P = op.raw[1] + α*op.raw[2] + α^2*op.raw[3];
                     q[op.keep] .= P[op.keep, op.keep] \ rhs; q)
function incident(N, alpha_f)
    opf = cavity_pencil(b, R, bf, om; N=N, pressure=:PN2)
    return solve(opf, alpha_f, stator_forcing(opf)[opf.keep])
end
function roughness(N, kappa)
    oph = cavity_pencil(b, R, bh, 0.0; N=N, pressure=:PN2)
    return roughness_response(oph, kappa)   # full field incl. wall values (needed for D q_h in the source)
end
function normalise(opt, x, zstar)
    f = full_state(opt, x); x ./ bary_interp(opt.z, f[opt.n+1:2opt.n], zstar)
end

# fixed normalisation height
opt, an, x, _ = setup(199); f = full_state(opt, x); nn = opt.n
zz = range(0, 0.3, length=3001)
zstar = zz[argmax([abs(bary_interp(opt.z, f[nn+1:2nn], t)) for t in zz])]
say(@sprintf("目标模态 α_n = %.8f %+.8fi（N = 199）；x 按 v(z* = %.4f) = 1 归一（z* 为 |v| 最大处，固定不变）。\n\n", real(an), imag(an), zstar))

# ------------------------------------------------------------------ V2
say("## V2 复数失谐检验（α_f = 0.1，N = 149）\n\n")
let N = 149, af = 0.1
    opt, an, x, y = setup(N); x = normalise(opt, x, zstar)
    qf = incident(N, af)
    rhs(αt) = (κ = αt - af; qh = roughness(N, κ);
               scattering_rhs(opt, αt, qh, κ, bh, qf, af, bf).total)
    Λ = residue_coefficient(opt, an, x, y, rhs(an))
    say(@sprintf("Λ_hf = %.6e %+.6ei\n\n| |δ| | ‖δ q_hf − x Λ‖/‖x Λ‖ |\n|---|---|\n", real(Λ), imag(Λ)))
    for d in (1e-2, 1e-3, 1e-4, 1e-5)
        δ = d*exp(im*π/4); αt = an + δ
        q = (opt.L0 + αt*opt.L1 + αt^2*opt.L2) \ rhs(αt)
        say(@sprintf("| %.0e | %.3e |\n", d, norm(δ .* q - x .* Λ)/norm(x .* Λ)))
    end
end
say("\n")

# ------------------------------------------------------------------ V3
say("## V3 网格收敛\n\n每单位静盘强迫幅值、单位粗糙元 Fourier 幅值（高度以间隙 h 为单位）。\n\n")
say("| α_f | N | |Λ_hf| | |Λ_vol| | |Λ_wall| | |Λ_dir| | |q_f|(转盘外缘) |\n|---|---|---|---|---|---|---|\n")
rows = String[]
for af in (0.05, 0.1, 0.2), N in (99, 149, 199, 249); let
    opt, an, x, y = setup(N); x = normalise(opt, x, zstar)
    qf = incident(N, af); κ = an - af; qh = roughness(N, κ)
    r = scattering_rhs(opt, an, qh, κ, bh, qf, af, bf)
    Λ, Λv, Λw = (residue_coefficient(opt, an, x, y, v) for v in (r.total, r.volume, r.wall))
    Λd = residue_coefficient(opt, an, x, y, stator_forcing(opt)[opt.keep])
    n = opt.n; edge = sqrt(sum(abs2(bary_interp(opt.z, qf[(j-1)*n+1:j*n], 0.085)) for j in 1:3))
    say(@sprintf("| %.2f | %d | %.6e | %.6e | %.6e | %.6e | %.6e |\n", af, N, abs(Λ), abs(Λv), abs(Λw), abs(Λd), edge))
    push!(rows, join((af, N, real(an), imag(an), real(Λ), imag(Λ), real(Λv), imag(Λv), real(Λw), imag(Λw), real(Λd), imag(Λd), edge), ","))
end; end
open(joinpath(OUT, "lambda_convergence.csv"), "w") do fh
    println(fh, "alpha_f,N,alpha_n_r,alpha_n_i,Lhf_r,Lhf_i,Lvol_r,Lvol_i,Lwall_r,Lwall_i,Ldir_r,Ldir_i,qf_edge")
    foreach(r -> println(fh, r), rows)
end
close(io)
println("\nwrote ", OUT)

#=
compute_cavity_scattering_scan.jl
=================================

alpha_f scan of the two excitation paths of a rotor traveling mode in the rotor–stator cavity
(CavityLST.jl, P_N-P_{N-2}, Re_h = 1000, a_s = 0, local parallel, rotor frame):

  direct    : stator-layer forcing at the target wavenumbers (alpha_n, n_t, omega_t)
              Lambda_dir = y^T f_s / Q
  scattering: stator forcing at (alpha_f, n_f, omega_t) x rotor-steady roughness (kappa, n_h),
              kappa = alpha_n - alpha_f (complex), n_h = n_t - n_f
              Lambda_hf = y^T b_hf(alpha_n) / Q   (volume cross source + wall Taylor term)
Both per unit forcing amplitude; roughness per unit Fourier amplitude with height in gap units.
Parity roughness height (flat forcing spectrum): h* = |Lambda_dir|/|Lambda_hf| (gap units),
h*/delta = h* sqrt(Re_h).

Targets : Type-I  R = 330, n_t = 16, omega_bar = 8   (alpha_n = 0.44878 - 0.01460i, growing)
          Type-II R = 100, n_t = -7, omega_bar = 8   (alpha_n = 0.27954 + 0.00062i, near neutral)
Incident: n_f = -8, -9, -10, -11  <->  c = -omega_bar/n_f = 1, 0.889, 0.8, 0.727 (all in [1-G_c, 1]);
          forcing direction e_r, e_theta, e_z (Gaussian at z_s = 0.955, sigma = 0.02);
          alpha_f in [-0.3, 0.42].
Contamination indicators per point:
  d_f = min |alpha_f - alpha_j| over eigenvalues of the incident pencil (beta_f, omega_t) -> q_f near-resonance
  d_h = min |kappa - kappa_j| over eigenvalues of the roughness pencil (beta_h, 0)      -> q_h self-resonance
  |q_f| at the rotor-layer edge, max |q_h|, and alpha_f R (local-parallel validity, want |alpha_f| R >~ 30).
Output: results/cavity_scattering_scan/scan.csv, report.md
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
include(joinpath(@__DIR__, "local_mode_tracking.jl"))
using LinearAlgebra, Printf
BLAS.set_num_threads(min(8, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "cavity_scattering_scan"); mkpath(OUT)
const Re_h, N = 1000.0, 149
const ZS, SIG, Z_EDGE = 0.955, 0.02, 0.085
const TARGETS = [(name="TypeI_R330_n16", R=330.0, n_t=16, wbar=8.0, alpha0=0.448780 - 0.014605im),
                 (name="TypeII_R100_n-7", R=100.0, n_t=-7, wbar=8.0, alpha0=0.279541 + 0.000615im)]
const NF = (-8, -9, -10, -11)
const AF = collect(-0.3:0.005:0.42)
b = cavity_baseflow(Re_h; N=300)

forcing(op, k) = (f = zeros(ComplexF64, 4op.n); f[(k-1)*op.n+1:k*op.n] .= exp.(-((op.z .- ZS).^2)./(2SIG^2)); f)
pmat(op, α) = op.raw[1] + α*op.raw[2] + α^2*op.raw[3]
function solve_full(op, α, rhs)
    q = zeros(ComplexF64, 4op.n); P = pmat(op, α); q[op.keep] .= P[op.keep, op.keep] \ rhs; q
end
velamp(op, q, z) = sqrt(sum(abs2(bary_interp(op.z, q[(j-1)*op.n+1:j*op.n], z)) for j in 1:3))
eigs_of(op) = [m.alpha for m in cavity_spectrum(op; shift=0.2+0im, maxabs=3.0)]

rows = NamedTuple[]
io = open(joinpath(OUT, "report.md"), "w")
say(s) = (print(s); print(io, s); flush(io))
say("# 腔体两条激发路径的 α_f 扫描\n\nRe_h = 1000，N = $(N)，P_N–P_{N−2}；单位静盘源强，粗糙元单位 Fourier 幅值（高度以间隙为单位）。\n\n")
for T in TARGETS
    R = T.R; om = T.wbar/R; bt = T.n_t/R
    opt = cavity_pencil(b, R, bt, om; N=N, pressure=:PN2)
    m = refine_mode(opt.L0, opt.L1, opt.L2, T.alpha0); an = m.alpha
    y = refine_mode(transpose(opt.L0), transpose(opt.L1), transpose(opt.L2), an).vector
    xf = full_state(opt, m.vector); n = opt.n
    zz = range(0, 0.3, length=3001); vz(t) = bary_interp(opt.z, xf[n+1:2n], t)
    zstar = zz[argmax(abs.(vz.(zz)))]; x = m.vector ./ vz(zstar)
    Λdir = [residue_coefficient(opt, an, x, y, forcing(opt, k)[opt.keep]) for k in 1:3]
    say(@sprintf("## %s：α_n = %.6f %+.6fi，x 按 v(z*=%.4f)=1；|Λ_dir| (e_r, e_θ, e_z) = %.3e, %.3e, %.3e\n\n",
                 T.name, real(an), imag(an), zstar, abs.(Λdir)...))
    for nf in NF
        nh = T.n_t - nf; bf = nf/R; bh = nh/R; c = -T.wbar/nf
        opf = cavity_pencil(b, R, bf, om; N=N, pressure=:PN2)
        oph = cavity_pencil(b, R, bh, 0.0; N=N, pressure=:PN2)
        ef = eigs_of(opf); eh = eigs_of(oph)
        t0 = time()
        for af in AF
            κ = an - af
            qh = roughness_response(oph, κ)
            dh, jh = findmin(abs.(eh .- κ)); df, jf = findmin(abs.(ef .- af))
            qhmax = maximum(abs.(qh[1:3n]))
            for k in 1:3
                qf = solve_full(opf, af, forcing(opf, k)[opf.keep])
                r = scattering_rhs(opt, an, qh, κ, bh, qf, af, bf)
                Λ, Λv, Λw = (residue_coefficient(opt, an, x, y, v) for v in (r.total, r.volume, r.wall))
                push!(rows, (target=T.name, R=R, n_t=T.n_t, wbar=T.wbar, alpha_n_r=real(an), alpha_n_i=imag(an),
                             n_f=nf, n_h=nh, c=c, dir=("r", "theta", "z")[k], alpha_f=af, kappa_r=real(κ), kappa_i=imag(κ),
                             absL_hf=abs(Λ), absL_vol=abs(Λv), absL_wall=abs(Λw), absL_dir=abs(Λdir[k]),
                             L_hf_r=real(Λ), L_hf_i=imag(Λ), hstar_gap=abs(Λdir[k])/abs(Λ),
                             hstar_delta=abs(Λdir[k])/abs(Λ)*sqrt(Re_h),
                             qf_edge=velamp(opf, qf, Z_EDGE), qh_max=qhmax,
                             d_f=df, near_f_r=real(ef[jf]), near_f_i=imag(ef[jf]),
                             d_h=dh, near_h_r=real(eh[jh]), near_h_i=imag(eh[jh]), alphaR=abs(af)*R))
            end
        end
        @printf("  %s n_f=%d (c=%.3f, n_h=%d): %.1f s\n", T.name, nf, c, nh, time() - t0)
    end
end
open(joinpath(OUT, "scan.csv"), "w") do f
    cols = collect(keys(rows[1])); println(f, join(string.(cols), ","))
    for r in rows; println(f, join([string(getproperty(r, c)) for c in cols], ",")); end
end
# summary: best clean point per target / n_f / direction (clean: d_f > 0.02, d_h > 0.02, |alpha_f| R >= 30)
say("## 每个 (目标, n_f, 方向) 的最小平价高度 h*/δ（干净点：d_f>0.02、d_h>0.02、|α_f|R≥30）\n\n")
say("| 目标 | n_f | c | 方向 | α_f | |Λ_hf| | |Λ_vol|/|Λ_wall| | |Λ_dir| | h*/δ | d_f | d_h |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
for T in TARGETS, nf in NF, d in ("r", "theta", "z")
    sel = [r for r in rows if r.target == T.name && r.n_f == nf && r.dir == d && r.d_f > 0.02 && r.d_h > 0.02 && r.alphaR >= 30]
    isempty(sel) && (say("| $(T.name) | $nf | | $d | 无干净点 | | | | | | |\n"); continue)
    r = sel[argmin([s.hstar_delta for s in sel])]
    say(@sprintf("| %s | %d | %.3f | %s | %.3f | %.3e | %.2f | %.3e | %.3e | %.3f | %.3f |\n", T.name, nf, r.c, d, r.alpha_f,
                 r.absL_hf, r.absL_vol/r.absL_wall, r.absL_dir, r.hstar_delta, r.d_f, r.d_h))
end
close(io)
println("wrote ", OUT)

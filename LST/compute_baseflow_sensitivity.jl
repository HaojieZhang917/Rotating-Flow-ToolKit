#=
compute_baseflow_sensitivity.jl
===============================

Step (b): sensitivity of the rotor target mode and of the two excitation paths to the base flow
(laminar Batchelor vs turbulent-stator models), Re_h = 1000, P_N-P_{N-2}, N = 149,
perturbation operator with molecular viscosity only (eddy viscosity only shapes the mean flow).

Targets (same as compute_cavity_scattering_scan.jl): Type-I R = 330, n_t = 16, omega_bar = 8;
Type-II R = 100, n_t = -7, omega_bar = 8. The turbulent base flow is evaluated at the target radius.
alpha_n is followed from the laminar mode along the homotopy  b(t) = (1-t) b_lam + t b_turb, t = 0 ... 1
(eigenvalue continuation; a full spectrum at t = 1 confirms the identification).
Incident n_f = -9, -10 (c = 0.889, 0.8), forcing e_r; alpha_f in [-0.3, 0.42].
Reports |Lambda_dir|, |Lambda_hf| and the parity height h*/delta = |Lambda_dir|/|Lambda_hf| sqrt(Re_h).
Output: results/baseflow_sensitivity/{report.md, scan.csv}
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
include(joinpath(@__DIR__, "local_mode_tracking.jl"))
using LinearAlgebra, Printf
BLAS.set_num_threads(min(8, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "baseflow_sensitivity"); mkpath(OUT)
const Re_h, N, NB = 1000.0, 149, 300
const ZS, SIG = 0.955, 0.02
const TARGETS = [(name="TypeI_R330_n16", R=330.0, n_t=16, wbar=8.0, alpha0=0.448780 - 0.014605im),
                 (name="TypeII_R100_n-7", R=100.0, n_t=-7, wbar=8.0, alpha0=0.279541 + 0.000615im)]
const NF = (-9, -10)
const AF = collect(-0.3:0.01:0.42)
const MODELS = (:laminar, :mixing_length, :sa)

io = open(joinpath(OUT, "report.md"), "w")
say(s) = (print(s); print(io, s); flush(io))
blend(a, b, t) = CavityBase(a.Re_h, a.a_s, a.z, (1-t).*a.F .+ t.*b.F, (1-t).*a.G .+ t.*b.G, (1-t).*a.H .+ t.*b.H,
                            (1-t)*a.Pi + t*b.Pi, a.D, 0.0, 0)
forcing(op) = (f = zeros(ComplexF64, 4op.n); f[1:op.n] .= exp.(-((op.z .- ZS).^2)./(2SIG^2)); f)
function solve_full(op, α, rhs)
    q = zeros(ComplexF64, 4op.n); P = op.raw[1] + α*op.raw[2] + α^2*op.raw[3]
    q[op.keep] .= P[op.keep, op.keep] \ rhs; q
end

lam = cavity_baseflow(Re_h; N=NB)
# turbulent base flows at the target radii; SA is continued in R from R = 100 (previous chi as seed)
bases = Dict{Tuple{Symbol,Float64},Any}()
for T in TARGETS
    bases[(:laminar, T.R)] = (base=lam, K=bary_interp(lam.z, lam.G, 0.5))
    bases[(:mixing_length, T.R)] = cavity_baseflow_turb(Re_h, T.R; N=NB)
end
let chi = nothing
    for R in vcat(collect(100.0:25.0:325.0), 330.0)
        r = cavity_baseflow_turb(Re_h, R; N=NB, model=:sa, chi0=chi); chi = r.chi
        any(T -> T.R == R, TARGETS) && (bases[(:sa, R)] = r)
    end
end

say("# 基流敏感性：层流 vs 湍流静盘层（混合长 / SA）\n\nRe_h = 1000，N = $(N)，P_N–P_{N−2}，扰动算子只含分子黏性；强迫 e_r；α_n 由层流模态沿同伦 b(t) 延拓。\n\n")
rows = String[]
for T in TARGETS
    R = T.R; om = T.wbar/R; bt = T.n_t/R
    say("## $(T.name)\n\n| 基流 | K | 转盘系核心速度 1−K | α_n | 延拓最大跳跃 |\n|---|---|---|---|---|\n")
    for M in MODELS
        bb = bases[(M, R)]; bturb = bb.base
        α = T.alpha0; jump = 0.0; x = nothing; opt = nothing
        for t in range(0, 1, length=21)
            bt_ = blend(lam, bturb, t)
            opt = cavity_pencil(bt_, R, bt, om; N=N, pressure=:PN2)
            m = refine_mode(opt.L0, opt.L1, opt.L2, α); jump = max(jump, abs(m.alpha - α)); α = m.alpha; x = m.vector
        end
        # identification check against the full spectrum at t = 1
        sp = cavity_spectrum(opt; shift=α); dmin = minimum(abs(s.alpha - α) for s in sp)
        y = refine_mode(transpose(opt.L0), transpose(opt.L1), transpose(opt.L2), α).vector
        xf = full_state(opt, x); n = opt.n
        zz = range(0, 0.3, length=3001); vz(t) = bary_interp(opt.z, xf[n+1:2n], t)
        zstar = zz[argmax(abs.(vz.(zz)))]; x = x ./ vz(zstar)
        Λdir = residue_coefficient(opt, α, x, y, forcing(opt)[opt.keep])
        say(@sprintf("| %s | %.4f | %.4f | %.6f %+.6fi | %.3f（全谱最近 %.1e） |\n", M, bb.K, 1 - bb.K, real(α), imag(α), jump, dmin))
        for nf in NF
            nh = T.n_t - nf; bf = nf/R; bh = nh/R
            opf = cavity_pencil(bturb, R, bf, om; N=N, pressure=:PN2)
            oph = cavity_pencil(bturb, R, bh, 0.0; N=N, pressure=:PN2)
            for af in AF
                κ = α - af
                qh = roughness_response(oph, κ)
                qf = solve_full(opf, af, forcing(opf)[opf.keep])
                Λ = residue_coefficient(opt, α, x, y, scattering_rhs(opt, α, qh, κ, bh, qf, af, bf).total)
                push!(rows, join((T.name, M, bb.K, real(α), imag(α), nf, -T.wbar/nf, af, abs(Λdir), abs(Λ),
                                  abs(Λdir)/abs(Λ)*sqrt(Re_h)), ","))
            end
        end
    end
    say("\n")
end
open(joinpath(OUT, "scan.csv"), "w") do f
    println(f, "target,model,K,alpha_n_r,alpha_n_i,n_f,c,alpha_f,absL_dir,absL_hf,hstar_delta"); foreach(r -> println(f, r), rows)
end
# summary table: h*/delta at selected alpha_f
say("## h*/δ 在若干 α_f 处的对比（n_f = −10，e_r）\n\n| 目标 | 基流 | |Λ_dir| | α_f=−0.2 | −0.1 | 0.1 | 0.2 | 0.3 |\n|---|---|---|---|---|---|---|---|\n")
for T in TARGETS, M in MODELS
    sel = [split(r, ",") for r in rows if startswith(r, T.name * "," * string(M) * ",")]
    g(af) = (k = findfirst(s -> parse(Int, s[6]) == -10 && abs(parse(Float64, s[8]) - af) < 1e-9, sel); k === nothing ? NaN : parse(Float64, sel[k][11]))
    ld = parse(Float64, sel[1][9])
    say(@sprintf("| %s | %s | %.3e | %.2e | %.2e | %.2e | %.2e | %.2e |\n", T.name, M, ld, g(-0.2), g(-0.1), g(0.1), g(0.2), g(0.3)))
end
close(io)
println("wrote ", OUT)

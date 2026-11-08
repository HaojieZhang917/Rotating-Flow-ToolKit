"""
compute_traveling_receptivity.jl
================================

Vortex–roughness receptivity coefficients of TRAVELING modes in the von Kármán
boundary layer (rotating frame, local parallel model). Roughness alone (omega = 0)
cannot excite these modes; the free-stream vortex supplies the frequency.

Euler constraint (inviscid frozen convection in the quiescent far field):
    omega_f = -beta_f  <=>  omega1_t = -m_f
Sum rules for the second-order response:
    m_t = m_f + n,   omega_t = omega_f,   alpha_t = alpha_f + alpha_h
so for a target (R, m_t, omega1_t):  m_f = -omega1_t,  n = m_t - m_f,  alpha_h = Re alpha_n - alpha_f.

Targets (from results/traveling_growth_maps):
    Type-II, low Re   : omega1 = 8, R = 200, m_t = -2
    slow Type-I       : omega1 = 4, R = 400, m_t = 22

Incident polarisations at beta_f != 0 (the two vortical far-field roots are degenerate):
    toroidal : w = 0                          (-> azimuthal swirl when beta_f -> 0)
    poloidal : alpha_f v - beta_f u = 0       (-> meridional when beta_f -> 0)
normalised by the velocity amplitude of the inviscid far-field solution at the
boundary-layer edge z_e = delta_99 (lab-frame G = 0.01).

Coefficient (same discretisation as `solve_coupled_response`: affine grid, Dirichlet edge):
    q_keep ≈ x * Lambda / (alpha_t - alpha_n),   Lambda = y^T b / Q,  Q = y^T (L1 + 2 alpha_n L2) x
x normalised v(z*) = 1 at the |v| peak of the target eigenfunction. Per complex
Fourier component; unit roughness amplitude.

Lambda is the exact residue: b(alpha) is continued analytically to the complex alpha_n,
i.e. the roughness distortion is solved at the complex kappa = alpha_n - alpha_f with the
quadratic pencil assembled directly (`qh_complex`). The residue is validated with complex
detuning alpha = alpha_n + delta on the same discrete operator (delta*q -> x*Lambda).
The direct, adjoint and forced problems share one discretisation (affine grid, Dirichlet
edge); truncation is checked by increasing zmax with N scaled proportionally.

The steady roughness distortion q_h(alpha_h, n/R) may itself be near a STATIONARY
mode at beta = n/R; ‖q_h‖ is recorded to flag that contamination.

Output: results/traveling_receptivity/
"""

include(joinpath(@__DIR__, "VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__, "solve_six_state_lns.jl"))
using .SixStateLNS
include(joinpath(@__DIR__, "local_mode_tracking.jl"))
using LinearAlgebra
using Printf
BLAS.set_num_threads(min(16, Sys.CPU_THREADS))

const OUT = joinpath(@__DIR__, "results", "traveling_receptivity"); mkpath(OUT)
const NMAIN = 129
const ZMAIN = 40.0
const ALPHA_F_LIST = vcat(collect(-0.80:0.05:-0.05), [-0.02, 0.02], collect(0.05:0.05:0.80))
const TARGETS = [
    (name="TypeII_w8", label="Type-II（低 Re）", omega1=8.0, R=200.0, m_t=-2, seed_omega1=7.9, seed=0.185 - 0.0161im),
    (name="TypeI_w4", label="慢速 Type-I", omega1=4.0, R=400.0, m_t=22, seed_omega1=4.0, seed=0.354 - 0.0365im),
]

base = solve_von_karman_baseflow()
const ZREF = let g = base.v .+ 1.0, j = findfirst(<(0.01), g)
    base.z[j-1] + (0.01 - g[j-1])*(base.z[j] - base.z[j-1])/(g[j] - g[j-1])
end
mkp(alpha, beta, omega; N, R) = LSTParameters(mode=:spatial, R=Float64(R), alpha=Float64(alpha), beta=Float64(beta),
                                              omega=Float64(omega), N=N, nev=8, spatial_shift=0.38+0im)

function pencil(R, beta, omega; N=NMAIN, zmax=ZMAIN, grid=:affine)
    p = mkp(0.38, beta, omega; N=N, R=R)
    b = assemble_spatial_blocks(base; p=p, grid=grid, zmax=zmax)
    L0, L1, L2 = VonKarmanLST._spatial_pencil(b, p)
    k = b.keep
    (L0=Matrix{ComplexF64}(L0[k, k]), L1=Matrix{ComplexF64}(L1[k, k]), L2=Matrix{ComplexF64}(L2[k, k]),
     raw=(Matrix{ComplexF64}(L0), Matrix{ComplexF64}(L1), Matrix{ComplexF64}(L2)), blocks=b)
end

"""Target mode on the forced-solver discretisation; seed continued in omega1 on the eigen grid first."""
function target_mode(t; N=NMAIN, zmax=ZMAIN, zstar=nothing, guess=nothing)
    R, beta = t.R, t.m_t/t.R
    a0 = if guess === nothing
        path = continue_mode(w -> pencil(R, beta, w/R; N=109, grid=:finite, zmax=30.0),
                             range(t.seed_omega1, t.omega1, length=6), t.seed)
        length(path) == 6 || error("seed continuation failed for $(t.name)")
        path[end].alpha
    else
        guess
    end
    ap = pencil(R, beta, t.omega1/R; N=N, zmax=zmax)
    dm = refine_mode(ap.L0, ap.L1, ap.L2, a0); am = refine_mode(transpose(ap.L0), transpose(ap.L1), transpose(ap.L2), dm.alpha)
    (dm.converged && am.converged && abs(dm.alpha - am.alpha) < 1e-10) || error("direct/adjoint eigenvalues disagree")
    n = N + 1; k = ap.blocks.keep; bc = collect(setdiff(1:4n, k))
    full = zeros(ComplexF64, 4n); full[k] .= dm.vector
    # closed-form Chebyshev–Gauss–Lobatto barycentric weights; the generic product formula overflows for N ≳ 200
    zb = Float64.(ap.blocks.z); bw = [(-1.0)^j*(j == 0 || j == N ? 0.5 : 1.0) for j in 0:N]; vb = full[n+1:2n]
    vint(zz) = VonKarmanLST._barycentric_interpolate(zb, vb, bw, zz)
    zstar === nothing && (zf = range(0.0, 15.0; length=6001); zstar = zf[argmax(abs.(vint.(zf)))])
    x = dm.vector ./ vint(zstar); y = am.vector
    Q = transpose(y)*(ap.L1 + 2dm.alpha*ap.L2)*x
    return (alpha=dm.alpha, x=x, y=y, Q=Q, L0=ap.raw[1], L1=ap.raw[2], L2=ap.raw[3], keep=k, bc=bc, n=n,
            N=N, zmax=zmax, zstar=zstar, R=R, beta=beta, omega=t.omega1/R, t=t,
            direct_residual=dm.residual, adjoint_residual=am.residual)
end

"""Incident vortical states (toroidal / poloidal) at (alpha_f, beta_f, omega_f = -beta_f)."""
function incident_states(alpha_f, beta_f, R; N=NMAIN, zmax=ZMAIN)
    p = mkp(alpha_f, beta_f, -beta_f; N=N, R=R)
    H = base.w[end]; k2 = alpha_f^2 + beta_f^2
    kz0 = -im*(H + sqrt(H^2 + 4k2))/2          # Doppler and Coriolis vanish: lambda^2 - H lambda - k^2 = 0
    near = sort(farfield_six_state_modes(base; p=p), by=m -> abs(m.kz - kz0))[1:2]
    maximum(abs(m.kz - kz0) for m in near) < 1e-8 || error("incident root mismatch at alpha_f=$(alpha_f)")
    S = hcat([m.state for m in near]...)
    pol = Dict("toroidal" => S*nullspace(reshape(S[3, :], 1, 2))[:, 1],
               "poloidal" => S*nullspace(reshape(alpha_f .* S[2, :] .- beta_f .* S[1, :], 1, 2))[:, 1])
    Dict(k => v ./ norm(v[1:3]) .* exp(im*kz0*(zmax - ZREF)) for (k, v) in pol), kz0, p
end

function first_order(tm, alpha_f, pol; N=tm.N, zmax=tm.zmax)
    mf = -tm.t.omega1
    st, kz0, p = incident_states(alpha_f, mf/tm.R, tm.R; N=N, zmax=zmax)
    as_four_state_response(solve_six_state_lns(base; p=p, zmax=zmax, boundary=:modal, incident_state=st[pol])), kz0
end

wall_rhs(tm, alpha, w) = (P = tm.L0 + alpha*tm.L1 + alpha^2*tm.L2; n = tm.n;
    qbc = zeros(ComplexF64, 4n); qbc[1] = w[1]; qbc[n+1] = w[2]; qbc[2n+1] = w[3]; -P[tm.keep, tm.bc]*qbc[tm.bc])

"""
Steady roughness distortion q_h at a possibly COMPLEX radial wavenumber kappa.
The four-variable operator is a quadratic pencil in alpha, so P_h(kappa) = L0 + kappa L1 + kappa^2 L2
is assembled directly (affine grid, wall q_h(0) = -[F'(0), G'(0), H'(0)], homogeneous edge),
exactly as `solve_roughness_distortion` does for real kappa.
"""
function qh_complex(R, beta_h, kappa; N, zmax)
    rp = pencil(R, beta_h, 0.0; N=N, zmax=zmax)
    L0, L1, L2 = rp.raw; n = N + 1; bc = collect(setdiff(1:4n, rp.blocks.keep))
    P = L0 + kappa*L1 + kappa^2*L2; rhs = zeros(ComplexF64, 4n)
    for row in bc
        P[row, :] .= 0; P[row, row] = 1
    end
    rhs[1] = -base.du[1]; rhs[n+1] = -base.dv[1]; rhs[2n+1] = -base.dw[1]
    p0 = mkp(0.0, beta_h, 0.0; N=N, R=R)
    q = P \ rhs
    return (q=q, D=rp.blocks.D, params=(R=p0.R, Ro=p0.Ro, Co=p0.Co, alpha=ComplexF64(kappa), beta=beta_h, omega=0.0),
            interior_residual=norm((L0 + kappa*L1 + kappa^2*L2)[rp.blocks.keep, :]*q, Inf))
end

"""Second-order unit forcing at the (possibly complex) target alpha_t = alpha_f + kappa."""
function forcing(tm, qf, kappa)
    nrough = tm.t.m_t + tm.t.omega1                 # n = m_t - m_f, m_f = -omega1
    qh = qh_complex(tm.R, nrough/tm.R, kappa; N=tm.N, zmax=tm.zmax)
    fc = nonlinear_cross_forcing(qh, qf, base; p_h=qh.params, p_f=qf.params)
    (abs(fc.target_beta - tm.beta) < 1e-13 && abs(fc.target_omega - tm.omega) < 1e-13) ||
        error("sum rule violated: beta $(fc.target_beta) vs $(tm.beta), omega $(fc.target_omega) vs $(tm.omega)")
    w = -wall_shear_response(qf)
    return (bvol=fc.source[tm.keep], bwall=wall_rhs(tm, fc.target_alpha, w), alpha_t=fc.target_alpha, qh=qh, wall=w)
end
"""
Far-field criterion rho = lambda_f / Re k_h. The incident vortex grows upward like e^{lambda_f z}
(lambda_f = -Im kz_inc); the irrotational far field of the steady distortion decays like
e^{-k_h z}, k_h = sqrt(kappa^2 + beta_h^2). For rho >= 1 the cross source does not decay above the
boundary layer and Lambda depends on the vertical extent of the domain (unbounded upper-deck interaction).
"""
farfield_rho(kz0, kappa, beta_h) = -imag(kz0)/real(sqrt(ComplexF64(kappa)^2 + beta_h^2))

coef(tm, f) = (vol = (transpose(tm.y)*f.bvol)/tm.Q; wall = (transpose(tm.y)*f.bwall)/tm.Q; (total=vol+wall, vol=vol, wall=wall))

report = IOBuffer()
say(s) = (println(s); println(report, s); flush(stdout))
const CONV_CONFIGS = [(249, 80.0), (129, 40.0), (169, 40.0), (129, 30.0), (189, 60.0)]   # first entry = reference
say("# 行进波的涡波–粗糙元感受性系数（von Kármán，局部平行，复数 α 留数）\n")
say(@sprintf("主网格 N=%d，zmax=%.0f（affine 网格、外缘齐次 Dirichlet；直接/伴随/受迫三者共用同一离散算子）。入射涡按无黏远场解在 z_e=δ99=%.4f 处单位速度归一；x 按 v(z*)=1 归一；单位粗糙度、每个复 Fourier 分量。Euler 约束 ω_f=-β_f，即 m_f=-ω̄_t；n=m_t-m_f。\n", NMAIN, ZMAIN, ZREF))
say("Λ 定义为严格留数：源项 b(α) 在复数 α_n 处解析延拓，即粗糙元畸变在复数 κ=α_n-α_f 处求解（四变量二次铅笔直接组装），Λ=yᵀb(α_n)/Q。表中 Λ_real 为旧做法（b 取在实数 α_t=Re α_n），用于说明修正量。\n")
allrows = NamedTuple[]
for t in TARGETS
    tm = target_mode(t)
    mf = -Int(t.omega1); nr = t.m_t - mf
    say("## $(t.label)：ω̄=$(t.omega1)，R=$(t.R)，m_t=$(t.m_t)（入射涡 m_f=$(mf)，粗糙元 n=$(nr)）\n")
    say(@sprintf("- α_n = `%.10f%+.10fi`（空间增长率 %.4e），β=%.6f，ω=%.6f，波角 %.2f°；直接/伴随残差 %.1e/%.1e；z*=%.3f。\n",
            real(tm.alpha), imag(tm.alpha), -imag(tm.alpha), tm.beta, tm.omega, atand(tm.beta, real(tm.alpha)),
            tm.direct_residual, tm.adjoint_residual, tm.zstar))
    # (a) complex-kappa distortion reproduces solve_roughness_distortion at real kappa
    kr = real(tm.alpha) - 0.2
    qa = qh_complex(tm.R, nr/tm.R, kr; N=tm.N, zmax=tm.zmax)
    qb = solve_roughness_distortion(base; p=mkp(0.0, 0.0, 0.0; N=tm.N, R=tm.R), kappa=kr, n=nr, target_beta=nr/tm.R,
                                    h_amplitude=1.0, zmax=tm.zmax, edge_state=zeros(ComplexF64, 4)).response
    say(@sprintf("- 一致性：实数 κ=%.4f 时 `qh_complex` 与 `solve_roughness_distortion` 相对差 %.1e；复数 κ 时内部残差 %.1e。\n",
            kr, norm(qa.q - qb.q)/norm(qb.q), qh_complex(tm.R, nr/tm.R, tm.alpha - 0.2; N=tm.N, zmax=tm.zmax).interior_residual))
    # (b) residue validation with complex detuning
    say("留数验证：α = α_n + δ（δ=|δ|e^{iπ/4}），同一离散算子直接求解 P(α)q = b(α)，检验 δ·q → x·Λ：\n")
    say("|偏振|α_f|\\|δ\\||\\|δ·q - x·Λ\\| / \\|x·Λ\\||\n|---|---:|---:|---:|")
    for pol in ("poloidal", "toroidal")
        af = 0.2; qf, _ = first_order(tm, af, pol)
        Λ = coef(tm, forcing(tm, qf, tm.alpha - af)).total
        for d in (1e-2, 1e-3, 1e-4)
            δ = d*exp(im*π/4); α = tm.alpha + δ
            f = forcing(tm, qf, α - af)
            P = (tm.L0 + α*tm.L1 + α^2*tm.L2)[tm.keep, tm.keep]
            q = P \ (f.bvol + f.bwall)
            say(@sprintf("|%s|%.2f|%.0e|%.2e|", pol, af, d, norm(δ .* q - tm.x .* Λ)/norm(tm.x .* Λ)))
        end
    end
    # (c) convergence in N and zmax
    say("\nN / zmax 收敛（参考 N=$(CONV_CONFIGS[1][1])、zmax=$(CONV_CONFIGS[1][2])；α_f=-0.2、0.4 满足 ρ<1，α_f=0.2 用于展示 ρ≥1 时的发散）：\n")
    say("|偏振|α_f|ρ|N|zmax|\\|Λ\\||arg Λ|相对差|\n|---|---:|---:|---:|---:|---:|---:|---:|")
    for pol in ("poloidal", "toroidal"), afc in (-0.2, 0.4, 0.2)
        ref = nothing
        for (N, L) in CONV_CONFIGS
            tt = (N, L) == (NMAIN, ZMAIN) ? tm : target_mode(t; N=N, zmax=L, zstar=tm.zstar, guess=tm.alpha)
            qfc, kzc = first_order(tt, afc, pol)
            Λ = coef(tt, forcing(tt, qfc, tt.alpha - afc)).total
            ref === nothing && (ref = Λ)
            say(@sprintf("|%s|%+.1f|%.2f|%d|%.0f|%.6e|%.5f|%.1e|", pol, afc, farfield_rho(kzc, tt.alpha - afc, nr/tt.R), N, L,
                    abs(Λ), angle(Λ), abs(Λ - ref)/abs(ref)))
        end
    end
    if false
    say("\nN / zmax 收敛（α_f=0.2；参考 N=$(CONV_CONFIGS[1][1])、zmax=$(CONV_CONFIGS[1][2])，N 随 zmax 同比增大以保持分辨率）：\n")
    say("|偏振|N|zmax|α_n|\\|Λ\\||arg Λ|相对差|\n|---|---:|---:|---|---:|---:|---:|")
    for pol in ("poloidal", "toroidal")
        ref = nothing
        for (N, L) in CONV_CONFIGS
            tt = (N, L) == (NMAIN, ZMAIN) ? tm : target_mode(t; N=N, zmax=L, zstar=tm.zstar, guess=tm.alpha)
            Λ = coef(tt, forcing(tt, first_order(tt, 0.2, pol)[1], tt.alpha - 0.2)).total
            ref === nothing && (ref = Λ)
            say(@sprintf("|%s|%d|%.0f|%.8f%+.8fi|%.6e|%.5f|%.1e|", pol, N, L, real(tt.alpha), imag(tt.alpha), abs(Λ), angle(Λ), abs(Λ - ref)/abs(ref)))
        end
    end
    end
    # (d) alpha_f scan
    for pol in ("poloidal", "toroidal"), af in ALPHA_F_LIST
        qf, kz0 = first_order(tm, af, pol)
        f = forcing(tm, qf, tm.alpha - af); c = coef(tm, f)
        cr = coef(tm, forcing(tm, qf, real(tm.alpha) - af)).total
        push!(allrows, (target=t.name, omega1=t.omega1, R=t.R, m_t=Float64(t.m_t), m_f=Float64(mf), n=Float64(nr),
                        alpha_n_re=real(tm.alpha), alpha_n_im=imag(tm.alpha), polarization=pol, alpha_f=af,
                        kappa_re=real(tm.alpha) - af, kappa_im=imag(tm.alpha), kz_inc_imag=imag(kz0),
                        Lambda_re=real(c.total), Lambda_im=imag(c.total), abs_Lambda=abs(c.total),
                        abs_vol=abs(c.vol), abs_wall=abs(c.wall), abs_Lambda_real=abs(cr),
                        real_vs_complex=abs(cr - c.total)/abs(c.total), qh_norm=norm(f.qh.q[1:3tm.n]),
                        rho_farfield=farfield_rho(kz0, tm.alpha - af, nr/tm.R)))
    end
    rs = [r for r in allrows if r.target == t.name]
    say("\nα_f 扫描（完整数据见 CSV）：\n")
    say("ρ=λ_f/Re k_h ≥ 1 的行标 *，该处上层相互作用不收敛，Λ 依赖计算域高度，不作为结论。\n")
    say("|偏振|α_f|Re κ|ρ|\\|Λ\\||\\|Λ_vol\\||\\|Λ_wall\\||\\|Λ_real\\||\\|Λ_real-Λ\\|/\\|Λ\\||‖q_h‖|\n|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for r in rs
        abs(r.alpha_f) in (0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8) || continue
        say(@sprintf("|%s|%+.2f%s|%+.4f|%.2f|%.4e|%.4e|%.4e|%.4e|%.2f|%.3f|", r.polarization, r.alpha_f, r.rho_farfield >= 1 ? "*" : "",
                r.kappa_re, r.rho_farfield, r.abs_Lambda, r.abs_vol, r.abs_wall, r.abs_Lambda_real, r.real_vs_complex, r.qh_norm))
    end
    for pol in ("poloidal", "toroidal")
        q = [r for r in rs if r.polarization == pol && abs(r.alpha_f) >= 0.1 && r.rho_farfield < 1]
        b = q[argmax([r.abs_Lambda for r in q])]
        say(@sprintf("\n- %s：|α_f|≥0.1 范围内 |Λ| 最大 %.4e，位于 α_f=%+.2f（Re κ=%+.4f）。", pol, b.abs_Lambda, b.alpha_f, b.kappa_re))
    end
    say("")
end
open(joinpath(OUT, "traveling_receptivity_scan.csv"), "w") do io
    cols = collect(keys(allrows[1])); println(io, join(string.(cols), ","))
    for r in allrows; println(io, join([string(getproperty(r, c)) for c in cols], ",")); end
end
say("## 解释边界\n")
say("- ρ=λ_f/Re k_h≥1 时，入射涡向上的增长快于粗糙元势流畸变的衰减，交叉源项在边界层外不衰减（对应 Wu 2001 三层结构中的上层贡献无界），局部模型中 Λ 不定义；物理上需要给出来流扰动在远场的有界剖面（例如定转子腔的有限间隙）。")
say("- 局部平行模型；Type-II 目标 λ_r/R≈0.17，非平行修正可能显著，结果作量级参考。")
say("- ‖q_h‖ 异常增大的 α_f 表示定常畸变本身接近 β=n/R 的定常模态（粗糙元一阶共振），该处 Λ 被一阶共振放大，不代表纯涡波散射。")
say("- 局部粗糙元 h(r)=(1/2π)∫ĥe^{iκr}dκ 激发的模态幅值为 i·ε_h·A_v·ĥ(α_n-α_f)·Λ，ĥ 在复数 κ 处取值；Λ 与定常 Type-I 的系数（results/receptivity_coefficient）使用同一归一化约定。")
write(joinpath(OUT, "traveling_receptivity_report.md"), String(take!(report)))
println("wrote ", OUT)

# Shared definitions for compute_lab_stationary_modes.jl, loaded on the master
# and on every Distributed worker. See that file's docstring for the method.

include(joinpath(@__DIR__, "VonKarmanLST.jl"))
using .VonKarmanLST
using LinearAlgebra
using Printf

BLAS.set_num_threads(1)

envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
envfloat(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d

const N_LIST = haskey(ENV, "N_LIST") ? parse.(Int, split(ENV["N_LIST"], ",")) : [49, 69, 89]
const ZMAX_SCAN = envfloat("ZMAX_SCAN", 30.0)
const ZMAX_CHECK = envfloat("ZMAX_CHECK", 40.0)
const TOL = envfloat("TOL", 1e-6)
const TOL_Z = envfloat("TOL_Z", 1e-4)
# Type-I 参考模态 dIm(alpha)/dsigma 约为 10，粗步长会跳到相邻分支，因此用几何细步长
const SIGMA_LIST = haskey(ENV, "SIGMA_LIST") ? parse.(Float64, split(ENV["SIGMA_LIST"], ",")) :
    [1e-3*1.6^k for k in 0:13]
const BB_JUMP = envfloat("BB_JUMP", 0.15)
const R_MIN = envfloat("R_MIN", 100.0)
const R_MAX = envfloat("R_MAX", 500.0)
const R_STEP = envfloat("R_STEP", 25.0)
const M_MIN = envint("M_MIN", 0)
const M_MAX = envint("M_MAX", 150)
const M_STEP = envint("M_STEP", 10)
const ALPHA_MAX = envfloat("ALPHA_MAX", 3.0)
const IMAG_MAX = envfloat("IMAG_MAX", 1.0)
const RES_MAX = envfloat("RES_MAX", 1e-8)

const OUTDIR = joinpath(@__DIR__, "results", get(ENV, "OUTNAME", "lab_stationary_modes"))
mkpath(OUTDIR)

base = solve_von_karman_baseflow()

# 实验室系周向速度 G = V + 1，以及 F、G 的二阶导（对 F'、G' 作有限差分），用于临界层与拐点
const ZB = base.z
const FB = base.u
const GB = base.v .+ 1.0
function _ddz(y, z)
    d = similar(y)
    d[1] = (y[2]-y[1])/(z[2]-z[1]); d[end] = (y[end]-y[end-1])/(z[end]-z[end-1])
    for j in 2:length(y)-1
        d[j] = (y[j+1]-y[j-1])/(z[j+1]-z[j-1])
    end
    d
end
const D2FB = _ddz(base.du, ZB)
const D2GB = _ddz(base.dv, ZB)

"""
旋转系 Doppler 剖面 ar*F + beta*(G-1) - omega 的第一个内部零点和最近拐点；不存在时返回 NaN。
默认 omega=-beta 时退化为实验室系剖面 ar*F + beta*G。
"""
function critical_layer(ar, beta, omega=-beta)
    u = ar .* FB .+ beta .* (GB .- 1.0) .- omega
    u2 = ar .* D2FB .+ beta .* D2GB
    zc = NaN
    for j in 1:length(u)-1
        if u[j]*u[j+1] < 0
            zc = ZB[j] - u[j]*(ZB[j+1]-ZB[j])/(u[j+1]-u[j]); break
        end
    end
    zi = Float64[]
    for j in 2:length(u2)-1
        u2[j]*u2[j+1] < 0 && push!(zi, ZB[j] - u2[j]*(ZB[j+1]-ZB[j])/(u2[j+1]-u2[j]))
    end
    zin = isempty(zi) ? NaN : (isnan(zc) ? zi[1] : zi[argmin(abs.(zi .- zc))])
    return zc, zin
end

"""sigma 为频率虚部：-i*(omega+i*sigma)*Ta = -i*omega*Ta + sigma*Ta。"""
function make_operator(R, m, omega; N, zmax=ZMAX_SCAN, sigma=0.0)
    p = LSTParameters(mode=:spatial, R=Float64(R), beta=Float64(m/R), omega=Float64(omega),
                      N=N, nev=4*(N+1), spatial_shift=0.38 + 0im)
    b = assemble_spatial_blocks(base; p=p, grid=:finite, zmax=zmax)
    L0, L1, L2 = VonKarmanLST._spatial_pencil(b, p)
    k = b.keep
    L0k = Matrix{ComplexF64}(L0[k,k])
    sigma == 0 || (L0k .+= sigma .* b.Ta[k,k])
    return (L0=L0k, L1=Matrix{ComplexF64}(L1[k,k]),
            L2=Matrix{ComplexF64}(L2[k,k]), z=b.z, keep=k, params=p)
end

function residual(op, a, q)
    r = (op.L0 + a*op.L1 + a^2*op.L2) * q
    return norm(r) / max(norm(op.L0)*norm(q), eps(Float64))
end

"""特征函数速度能量峰位置及 z>10 的能量占比（按节点求和的粗指标）。"""
function profile_metrics(op, q)
    n = length(op.z); full = zeros(ComplexF64, 4n); full[op.keep] .= q
    e = abs2.(full[1:n]) .+ abs2.(full[n+1:2n]) .+ abs2.(full[2n+1:3n])
    tot = sum(e)
    return op.z[argmax(e)], tot > 0 ? sum(e[op.z .> 10.0])/tot : NaN
end

function spectrum(R, m, omega, N; zmax=ZMAX_SCAN, sigma=0.0)
    op = make_operator(R, m, omega; N=N, zmax=zmax, sigma=sigma)
    modes = solve_spatial(op; shift=0.38+0im, nev=size(op.L0,1))
    out = NamedTuple[]
    for x in modes
        a = x.alpha
        (isfinite(a) && abs(real(a)) <= ALPHA_MAX && abs(imag(a)) <= IMAG_MAX) || continue
        # m=0 时 alpha=0 是压力规范零模（参见 coupled_validation 的 kappa=0 检查），不是物理模态
        abs(a) < 1e-4 && continue
        rr = residual(op, a, x.vector)
        rr < RES_MAX || continue
        push!(out, (alpha=a, vector=x.vector, residual=rr))
    end
    return op, out
end

reldiff(a, b) = abs(a-b)/max(abs(a), 0.05)
nearest(a, s) = isempty(s) ? (Inf, a) : (d = [reldiff(a, y.alpha) for y in s]; j = argmin(d); (d[j], s[j].alpha))

"""
沿 SIGMA_LIST 追踪 alpha（最低 N），返回 (branch, 最后可靠追踪点的 Im alpha, 该点 sigma)。

每步取最近特征值；若绝对跳跃超过 BB_JUMP，或次近特征值距离不到最近者的
2 倍，则视为歧义并停止。追踪到 sigma_max 时按 Im alpha 符号判定
（sigma_max 远大于本问题的时间增长率）；提前停止时要求 Im alpha 的符号与
其变化方向一致，否则记为 undetermined。
"""
function briggs_bers(alpha0, sigma_specs)
    a = alpha0; s_last = 0.0
    for (s, spec) in zip(SIGMA_LIST, sigma_specs)
        isempty(spec) && break
        d = sort([abs(y.alpha - a) for y in spec])
        (d[1] > BB_JUMP || (length(d) > 1 && d[2] < 2d[1])) && break
        a = spec[argmin([abs(y.alpha - a) for y in spec])].alpha
        s_last = s
    end
    dai = imag(a) - imag(alpha0)
    br = if s_last == SIGMA_LIST[end]
        imag(a) > 0 ? "plus" : "minus"
    elseif imag(a) > 0 && dai > 0
        "plus"
    elseif imag(a) < 0 && dai < 0
        "minus"
    else
        "undetermined"
    end
    return (br, imag(a), s_last)
end

"""以最高 N 为参考，返回通过 N 收敛和基流截断高度检验的模态列表。"""
function converged_modes(R, m, omega; branch=true)
    Ns = sort(N_LIST)
    specs = [spectrum(R, m, omega, N)[2] for N in Ns[1:end-1]]
    opref, ref = spectrum(R, m, omega, Ns[end])
    cand = [x for x in ref if nearest(x.alpha, specs[end])[1] < TOL]
    isempty(cand) && return NamedTuple[]
    zspec = spectrum(R, m, omega, Ns[end-1]; zmax=ZMAX_CHECK)[2]
    sspecs = branch ? [spectrum(R, m, omega, Ns[1]; sigma=s)[2] for s in SIGMA_LIST] : []
    out = NamedTuple[]
    for x in cand
        d = [nearest(x.alpha, s)[1] for s in specs]
        dz = nearest(x.alpha, zspec)[1]
        a_low = nearest(x.alpha, specs[1])[2]
        br, ai_s, s_last = branch ? briggs_bers(a_low, sspecs) : ("none", NaN, NaN)
        zc, zi = critical_layer(real(x.alpha), m/R, omega)
        zp, ff = profile_metrics(opref, x.vector)
        push!(out, (R=Float64(R), m=Float64(m), beta=m/R, omega=Float64(omega),
                    alpha_real=real(x.alpha), alpha_imag=imag(x.alpha), growth=-imag(x.alpha),
                    residual=x.residual, dN_mid=d[end], dN_low=d[1], dzmax=dz,
                    zmax_ok=dz < TOL_Z, branch=br, imag_alpha_bb=ai_s, sigma_bb=s_last,
                    beta_over_alpha=(m/R)/real(x.alpha), z_crit=zc, z_inflect=zi,
                    z_peak=zp, farfield_fraction=ff))
    end
    return out
end

function write_csv(path, rows)
    isempty(rows) && (write(path, ""); return)
    cols = collect(keys(rows[1]))
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for r in rows
            println(io, join([string(getproperty(r, c)) for c in cols], ","))
        end
    end
end


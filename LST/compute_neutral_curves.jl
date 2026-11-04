"""
compute_neutral_curves.jl
========================

在已有 `VonKarmanLST.jl` 空间二次特征值求解器上，计算固定全局频率
omega1 下的局部空间中性点。这里使用

    beta  = m/R,       omega = omega1/R,

并严格把中性条件定义为 Im(alpha)=0。程序先用低分辨率粗搜索寻找
候选分支，再在 R 方向作一维局部根求解，并沿 m 方向用 alpha 和速度
特征函数重叠进行连续追踪。它不会用“最小 Im(alpha)”自动冒充物理
Type-I/Type-II 模态；只有与已验证 Type-I 参考模态连续相连的分支才
自动标记为 Type-I，其余分支保留为 unclassified，等待物理参考模态。

高分辨率验证可通过命令行环境变量设置：
  NK_SCAN, N_SCAN, ZMAX_SCAN, M_MIN, M_MAX, M_STEP, R_MIN, R_MAX, R_STEP
"""

include(joinpath(@__DIR__, "VonKarmanLST.jl"))
using .VonKarmanLST
using LinearAlgebra
using Serialization

const OUTDIR = joinpath(@__DIR__, "results", "neutral_curves")
mkpath(OUTDIR)

envint(k, d) = haskey(ENV, k) ? parse(Int, ENV[k]) : d
envfloat(k, d) = haskey(ENV, k) ? parse(Float64, ENV[k]) : d

const N_SCAN = envint("N_SCAN", 29)
const ZMAX_SCAN = envfloat("ZMAX_SCAN", 30.0)
const M_MIN = envint("M_MIN", -80)
const M_MAX = envint("M_MAX", 150)
const M_STEP = envint("M_STEP", 10)
const R_MIN = envfloat("R_MIN", 50.0)
const R_MAX = envfloat("R_MAX", 800.0)
const R_STEP = envfloat("R_STEP", 50.0)
const ROOT_ITERS = envint("ROOT_ITERS", 12)

const FREQUENCIES = haskey(ENV, "OMEGA1_LIST") ?
    parse.(Float64, split(ENV["OMEGA1_LIST"], ",")) : [-5.0, 0.0, 4.0, 7.9, 10.0]

base = solve_von_karman_baseflow()

"""把局部 beta、omega 按全局 m、omega1 生成空间算子。"""
function make_operator(R, m, omega1; N=N_SCAN, zmax=ZMAX_SCAN)
    p = LSTParameters(mode=:spatial, R=Float64(R),
                      beta=Float64(m/R), omega=Float64(omega1/R),
                      N=N, nev=4*(N+1), spatial_shift=0.38 + 0im)
    b = assemble_spatial_blocks(base; p=p, grid=:finite, zmax=zmax)
    L0, L1, L2 = VonKarmanLST._spatial_pencil(b, p)
    k = b.keep
    return (L0=Matrix{ComplexF64}(L0[k,k]),
            L1=Matrix{ComplexF64}(L1[k,k]),
            L2=Matrix{ComplexF64}(L2[k,k]),
            z=b.z, D=b.D, params=p, blocks=b)
end

function residual(op, a, q)
    r = (op.L0 + a*op.L1 + a^2*op.L2) * q
    return norm(r) / max(norm(op.L0)*norm(q), eps(Float64))
end

function velocity_overlap(a, b)
    na, nb = norm(a), norm(b)
    na == 0 || nb == 0 ? 0.0 : abs(dot(conj(a), b))/(na*nb)
end

"""保留有限、低残差、主导径向波数在物理范围内的候选。"""
function physical_modes(op)
    all = solve_spatial(op; shift=0.38+0im, nev=size(op.L0,1))
    out = NamedTuple[]
    for x in all
        a = x.alpha
        (isfinite(real(a)) && isfinite(imag(a)) &&
         0.02 <= real(a) <= 1.5 && abs(imag(a)) <= 1.5) || continue
        rr = residual(op, a, x.vector)
        rr < 5e-7 || continue
        q = ComplexF64.(x.vector)
        q ./= max(norm(q), eps(Float64))
        push!(out, (alpha=a, vector=q, residual=rr))
    end
    return out
end

function select_mode(modes; alpha_ref=nothing, vector_ref=nothing)
    isempty(modes) && return nothing
    if alpha_ref === nothing
        return modes[argmin([abs(imag(x.alpha)) + 0.02*abs(real(x.alpha)-0.38)
                             for x in modes])]
    end
    scores = Float64[]
    for x in modes
        da = abs(x.alpha-alpha_ref)/(0.10+abs(alpha_ref))
        ov = vector_ref === nothing ? 0.0 : velocity_overlap(vector_ref, x.vector)
        push!(scores, da + 0.35*(1.0-ov))
    end
    return modes[argmin(scores)]
end

"""给定参考 alpha/特征函数，在一个参数点追踪同一离散分支。"""
function branch_mode(R, m, omega1; alpha_ref=nothing, vector_ref=nothing)
    op = make_operator(R, m, omega1)
    modes = physical_modes(op)
    x = select_mode(modes; alpha_ref=alpha_ref, vector_ref=vector_ref)
    x === nothing && return nothing
    return merge(x, (R=Float64(R), m=Float64(m), omega1=Float64(omega1),
                     beta=Float64(m/R), omega=Float64(omega1/R), op=op))
end

"""在粗 R 网格上寻找 Im(alpha) 的变号区间。"""
function bracket_for_m(m, omega1, previous)
    # 有前一个中性点时，从该点沿 R 小步继续追踪；这样不会在
    # α 谱交叉处因“最近 alpha”跳到另一条离散分支。
    if previous !== nothing
        x0 = branch_mode(previous.R, m, omega1;
                         alpha_ref=previous.alpha, vector_ref=previous.vector)
        x0 === nothing && return nothing
        vals = [x0]
        for direction in (-1, 1)
            last = x0
            R = previous.R + direction*min(R_STEP, 10.0)
            while R_MIN <= R <= R_MAX
                x = branch_mode(R, m, omega1;
                                alpha_ref=last.alpha, vector_ref=last.vector)
                x === nothing && break
                push!(vals, x)
                # 找到变号后无需继续走到另一条分支。
                if imag(x.alpha)*imag(last.alpha) <= 0
                    break
                end
                last = x
                R += direction*min(R_STEP, 10.0)
            end
        end
        sort!(vals, by=x->x.R)
    else
        Rs = collect(R_MIN:R_STEP:R_MAX)
        vals = NamedTuple[]
        for R in Rs
            x = branch_mode(R, m, omega1; alpha_ref=nothing, vector_ref=nothing)
            x === nothing && continue
            push!(vals, x)
        end
    end
    length(vals) < 2 && return nothing
    brackets = NamedTuple[]
    for j in 1:length(vals)-1
        a, b = vals[j], vals[j+1]
        if imag(a.alpha) == 0 || imag(a.alpha)*imag(b.alpha) <= 0
            score = previous === nothing ? abs(imag(a.alpha))+abs(imag(b.alpha)) :
                    abs((a.R+b.R)/2-previous.R)
            push!(brackets, (a=a, b=b, score=score))
        end
    end
    isempty(brackets) && return nothing
    return brackets[argmin([q.score for q in brackets])]
end

"""对一个 m 求中性根；每个二分点都使用相邻端点的模态信息。"""
function neutral_at_m(m, omega1, previous)
    br = bracket_for_m(m, omega1, previous)
    br === nothing && return nothing
    lo, hi = br.a, br.b
    for _ in 1:ROOT_ITERS
        Rmid = 0.5*(lo.R+hi.R)
        t = (Rmid-lo.R)/max(hi.R-lo.R, eps())
        aref = (1-t)*lo.alpha + t*hi.alpha
        vref = (1-t)*lo.vector + t*hi.vector
        mid = branch_mode(Rmid, m, omega1; alpha_ref=aref, vector_ref=vref)
        mid === nothing && break
        if imag(lo.alpha)*imag(mid.alpha) <= 0
            hi = mid
        else
            lo = mid
        end
    end
    x = abs(imag(lo.alpha)) < abs(imag(hi.alpha)) ? lo : hi
    err = abs(imag(x.alpha))
    # 根求解必须真正达到中性；离散残差很小不能替代 Im(alpha)=0。
    err < 2e-4 || return nothing
    return merge(x, (neutral_error=err,))
end

"""从 omega1=0 的已知 Malik Type-I 点构造参考模态。"""
function type1_reference()
    R0 = 285.36; m0 = R0*0.07759
    op = make_operator(R0, m0, 0.0)
    modes = physical_modes(op)
    x = modes[argmin([abs(a.alpha-(0.38505249+0.00034867im)) for a in modes])]
    return merge(x, (R=R0, m=m0, omega1=0.0, beta=m0/R0, omega=0.0, op=op))
end

"""对一个频率沿 m 正负方向追踪连续中性支。"""
function track_type1(omega1; seed_ref=nothing)
    ref = seed_ref === nothing ? type1_reference() : seed_ref
    m0 = ref.m
    # 先计算最接近参考点的整数/连续 m，再从该点分别向两侧延拓。
    # 已知参考点的 m 不应先四舍五入，否则可能落到另一条近邻谱分支。
    mseed = m0
    seed = neutral_at_m(mseed, omega1, ref)
    seed === nothing && return NamedTuple[]
    rows = NamedTuple[]
    push!(rows, merge(seed, (classification="Type I", branch="Type-I", tracking_overlap=velocity_overlap(ref.vector, seed.vector))))
    for direction in (-1, 1)
        prev = seed
        m = mseed + direction*M_STEP
        while M_MIN <= m <= M_MAX
            x = neutral_at_m(m, omega1, prev)
            x === nothing && break
            ov = velocity_overlap(prev.vector, x.vector)
            ov > 0.35 || break
            push!(rows, merge(x, (classification="Type I", branch="Type-I", tracking_overlap=ov)))
            prev = x
            m += direction*M_STEP
        end
    end
    return rows
end

"""在固定 m 的中性点上对 omega1 作小步连续延拓。"""
function continue_frequency(seed, omega_from, omega_to; step=0.5)
    seed === nothing && return nothing
    delta = omega_to - omega_from
    n = max(1, ceil(Int, abs(delta)/step))
    current = seed
    for j in 1:n
        w = omega_from + delta*j/n
        x = neutral_at_m(current.m, w, current)
        x === nothing && return nothing
        current = x
    end
    return current
end

function write_csv(path, rows)
    cols = [:frequency, :R, :m, :beta, :omega, :alpha_real, :alpha_imag,
            :neutral_error, :residual, :tracking_overlap, :classification, :branch]
    open(path, "w") do io
        println(io, join(string.(cols), ","))
        for r in rows
            vals = [r.frequency, r.R, r.m, r.beta, r.omega, real(r.alpha),
                    imag(r.alpha), r.neutral_error, r.residual,
                    r.tracking_overlap, r.classification, r.branch]
            println(io, join(string.(vals), ","))
        end
    end
end

allrows = NamedTuple[]
allmodes = Dict{String,Any}()
ref0 = type1_reference()
# 先计算零频率参考曲线；非零频率从该点以小步长延拓，避免
# 直接在相隔很大的频率之间跳到另一条空间谱分支。
computed = Dict{Float64,Vector{Any}}()
seed_by_frequency = Dict{Float64,Any}(0.0 => ref0)
freq_order = sort(unique(vcat([0.0], FREQUENCIES)), by=abs)
for w1 in freq_order
    if w1 != 0.0
        if w1 > 0
            prevw = maximum([w for w in keys(seed_by_frequency) if w >= 0 && w < w1]; init=0.0)
        else
            prevw = minimum([w for w in keys(seed_by_frequency) if w <= 0 && w > w1]; init=0.0)
        end
        seed_by_frequency[w1] = continue_frequency(seed_by_frequency[prevw], prevw, w1)
    end
    println("tracking omega1=", w1)
    seed = seed_by_frequency[w1]
    rows = seed === nothing ? NamedTuple[] : track_type1(w1; seed_ref=seed)
    for r in rows
        push!(allrows, merge(r, (frequency=w1, alpha=r.alpha)))
        key = "omega1=$(w1)_m=$(r.m)"
        allmodes[key] = (frequency=w1, R=r.R, m=r.m, beta=r.beta,
                         omega=r.omega, alpha=r.alpha, vector=r.vector,
                         z=r.op.z, N=N_SCAN, zmax=ZMAX_SCAN,
                         residual=r.residual, neutral_error=r.neutral_error,
                         classification=r.classification)
    end
    if !isempty(rows)
        m0idx = argmin([abs(r.m-ref0.m) for r in rows])
        seed_by_frequency[w1] = rows[m0idx]
    end
    println("  points=", length(rows))
end

write_csv(joinpath(OUTDIR, "neutral_type1_curves.csv"), allrows)
serialize(joinpath(OUTDIR, "neutral_direct_modes.jls"), allmodes)

open(joinpath(OUTDIR, "neutral_solver_audit.md"), "w") do io
    println(io, "# 空间 LST 中性曲线求解器核查")
    println(io, "")
    println(io, "- 算子：`VonKarmanLST.jl` 的 `L0 + alpha*L1 + alpha^2*L2 = 0`，未重新构造 LST。")
    println(io, "- 参考系：随旋转圆盘坐标系；基本流为 `U=F, V=G-1, W=H`。")
    println(io, "- 全局到局部转换：`beta=m/R`、`omega=omega1/R`。算子中的 `beta*R` 与 `omega` 项均按现有工程原样使用；没有再次乘 R。")
    println(io, "- 空间中性条件：`imag(alpha)=0`；增长率使用 `sigma=-imag(alpha)`，没有混用 temporal 特征值。")
    println(io, "- 边界条件：有限 `zmax=$(ZMAX_SCAN)` 的齐次壁面/外缘截断，采用当前已验证的谱离散入口。")
    println(io, "- 分支追踪：R 粗网格变号 + 一维二分；相邻点同时使用 alpha 连续性和离散速度特征函数重叠。")
    println(io, "- Type-I 参考：R=285.36、beta=0.07759、alpha约为 `0.38505249+0.00034867im` 的现有回归点。")
    println(io, "- Type-II 自动识别：当前没有将未有物理参考种子的分支强行命名为 Type-II；若没有独立 Type-II 参考模态，结果会标为未分类。")
    println(io, "")
    println(io, "## 本次配置")
    println(io, "`N=$(N_SCAN)`, `zmax=$(ZMAX_SCAN)`, `R=$(R_MIN):$(R_STEP):$(R_MAX)`, `m=$(M_MIN):$(M_STEP):$(M_MAX)`, `omega1=$(FREQUENCIES)`。")
    println(io, "")
    println(io, "输出数据见 `neutral_type1_curves.csv`，中性点直接特征函数见 `neutral_direct_modes.jls`。")
end

open(joinpath(OUTDIR, "neutral_curves_report.md"), "w") do io
    println(io, "# 不同全局频率下空间中性曲线计算报告")
    println(io, "")
    println(io, "本报告仅报告当前离散原始变量空间 LST 能够稳定追踪的分支。")
    println(io, "")
    for w1 in FREQUENCIES
        rr = [r for r in allrows if r.frequency == w1]
        println(io, "- `omega1=$(w1)`: $(length(rr)) 个 Type-I 候选中性点。")
        if !isempty(rr)
            println(io, "  - R 范围 `$(minimum(r.R for r in rr))` 到 `$(maximum(r.R for r in rr))`；m 范围 `$(minimum(r.m for r in rr))` 到 `$(maximum(r.m for r in rr))`。")
            println(io, "  - 最大 `|Im(alpha)|=$(maximum(r.neutral_error for r in rr))`，最大离散残差 `$(maximum(r.residual for r in rr))`。")
        end
    end
    println(io, "")
    println(io, "## 解释边界")
    println(io, "当前自动流程验证了 Type-I 参考分支的连续延拓；Type-II 需要独立的已确认直接模态种子，不能由最小 `Im(alpha)` 规则推断。")
    println(io, "频率非零时若没有连续分支，表示在给定 R、m 搜索范围或当前有限外缘算子下未找到可靠 Type-I 延拓，不等于物理上不存在模态。")
end

println("wrote neutral curves to ", OUTDIR)

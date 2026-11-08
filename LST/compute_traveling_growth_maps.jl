"""
compute_traveling_growth_maps.jl
================================

Growth-rate maps and neutral curves of traveling modes in the (R, m) plane at the
frequencies of Thomas (2024, Phys. Fluids 36, 014117, Fig. 2):
omega1 = omega*R in {0, -5, 4, 7.9, 10}, rotating frame, beta = m/R.

1. map: on a dense (R, m) grid every point is processed by `map_point`
   (traveling_map_core.jl): N_FIND full QZ for discovery, local refinement on
   N_MID/N_REF for convergence, local zmax check and local Briggs–Bers test.
   The map value is -Im(alpha) of the least-damped physical (plus-branch, converged,
   zmax-insensitive, far-field fraction < 0.1) mode.
2. neutral points: sign changes of Im(alpha) of the least-damped physical mode
   between neighbouring grid points with continuous alpha (|Δα| < 0.1), along R at
   fixed m and along m at fixed R (m is a continuous local parameter), are refined by
   bisection with the local solver on N_REF.

Only Re(alpha) > 0 is kept (POSITIVE_ALPHA=1, the convention of Thomas 2024); the
Re(alpha) < 0 modes are the images (m, alpha, omega) -> (-m, -conj(alpha), -omega) of the
same physical waves at -omega. All m are therefore scanned for every omega1.

Environment: NPROCS, OMEGA1_LIST, R_MIN, R_MAX, R_STEP, M_MIN, M_MAX, M_STEP, OUTNAME
Output: results/traveling_growth_maps/
"""

using Distributed
using LinearAlgebra
using Printf

haskey(ENV, "OUTNAME") || (ENV["OUTNAME"] = "traveling_growth_maps")
haskey(ENV, "R_MIN") || (ENV["R_MIN"] = "50")
haskey(ENV, "R_MAX") || (ENV["R_MAX"] = "800")
haskey(ENV, "R_STEP") || (ENV["R_STEP"] = "25")
haskey(ENV, "M_MIN") || (ENV["M_MIN"] = "-14")
haskey(ENV, "M_MAX") || (ENV["M_MAX"] = "80")
haskey(ENV, "M_STEP") || (ENV["M_STEP"] = "2")
const NPROCS = haskey(ENV, "NPROCS") ? parse(Int, ENV["NPROCS"]) : 16
nprocs() == 1 && NPROCS > 1 && addprocs(NPROCS)
@everywhere include($(joinpath(@__DIR__, "traveling_map_core.jl")))

const OMEGA1_LIST = haskey(ENV, "OMEGA1_LIST") ? parse.(Float64, split(ENV["OMEGA1_LIST"], ",")) :
    [0.0, -5.0, 4.0, 7.9, 10.0]

# ---- regression: omega=0 Type-I reference through the fast filter ----
BLAS.set_num_threads(min(16, Sys.CPU_THREADS))
const ALPHA_TYPE1_REF = 0.38505249 + 0.00034867im
reg = map_point(0.0, 285.36, 285.36*0.07759)
ireg = isempty(reg) ? 0 : argmin([abs(complex(r.alpha_real, r.alpha_imag) - ALPHA_TYPE1_REF) for r in reg])
reg_err = ireg == 0 ? Inf : abs(complex(reg[ireg].alpha_real, reg[ireg].alpha_imag) - ALPHA_TYPE1_REF)
reg_err < 1e-6 && map_physical(reg[ireg]) || error("regression failed: err=$(reg_err), row=$(ireg == 0 ? nothing : reg[ireg])")
println(@sprintf("regression ok: |dalpha|=%.2e, branch=%s", reg_err, reg[ireg].branch)); flush(stdout)
BLAS.set_num_threads(1)

Rs = collect(R_MIN:R_STEP:R_MAX)
ms = collect(M_MIN:M_STEP:M_MAX)
pts = [(w, R, m) for w in OMEGA1_LIST for R in Rs for m in ms]
t0 = time()
prog = RemoteChannel(() -> Channel{Int}(length(pts)))
printer = @async for k in 1:length(pts)
    take!(prog)
    if k % 200 == 0 || k == length(pts)
        el = time() - t0
        println(@sprintf("  map %d/%d points (%.0f%%), %.0f s elapsed, ETA %.0f s", k, length(pts),
                         100k/length(pts), el, el*(length(pts)-k)/k))
        flush(stdout)
    end
end
results = pmap(pts) do (w, R, m)
    r = map_point(w, R, m)
    put!(prog, 1)
    r
end
wait(printer)
allrows = reduce(vcat, results; init=NamedTuple[])
write_csv(joinpath(OUTDIR, "map_modes.csv"), allrows)

best = Dict{Tuple{Float64,Float64,Float64},NamedTuple}()
for (i, p) in enumerate(pts)
    rs = filter(map_physical, results[i])
    isempty(rs) || (best[Float64.(p)] = rs[argmin([r.alpha_imag for r in rs])])
end
write_csv(joinpath(OUTDIR, "map_least_damped.csv"), [best[Float64.(p)] for p in pts if haskey(best, Float64.(p))])
println("map done in $(round(time()-t0)) s; ", length(best), "/", length(pts), " points with a physical mode"); flush(stdout)

# ---- neutral points: brackets along R (fixed m) and along m (fixed R) ----
brackets = []
for w in OMEGA1_LIST
    for m in ms, j in 1:length(Rs)-1
        a = get(best, (w, Rs[j], Float64(m)), nothing); b = get(best, (w, Rs[j+1], Float64(m)), nothing)
        (a === nothing || b === nothing) && continue
        a.alpha_imag*b.alpha_imag < 0 && abs(complex(a.alpha_real - b.alpha_real, a.alpha_imag - b.alpha_imag)) < 0.1 &&
            push!(brackets, (w, :R, Float64(m), Rs[j], Rs[j+1], complex(a.alpha_real, a.alpha_imag)))
    end
    for R in Rs, j in 1:length(ms)-1
        a = get(best, (w, R, Float64(ms[j])), nothing); b = get(best, (w, R, Float64(ms[j+1])), nothing)
        (a === nothing || b === nothing) && continue
        a.alpha_imag*b.alpha_imag < 0 && abs(complex(a.alpha_real - b.alpha_real, a.alpha_imag - b.alpha_imag)) < 0.1 &&
            push!(brackets, (w, :m, R, Float64(ms[j]), Float64(ms[j+1]), complex(a.alpha_real, a.alpha_imag)))
    end
end
println("refining $(length(brackets)) neutral brackets ..."); flush(stdout)

@everywhere function refine_bracket(w, dir, fixed, p0, p1, a0)
    RM(p) = dir == :R ? (p, fixed) : (fixed, p)
    build(p) = (Rm = RM(p); make_operator(Rm[1], Rm[2], w/Rm[1]; N=N_REF))
    path = continue_mode(build, range(p0, p1, length=9), a0; max_jump=0.05)
    length(path) == 9 || return nothing
    j = findfirst(k -> imag(path[k].alpha)*imag(path[k+1].alpha) <= 0, 1:8)
    j === nothing && return nothing
    lo, hi = path[j], path[j+1]
    for _ in 1:30
        pm = (lo.param + hi.param)/2
        o = build(pm); x = refine_mode(o.L0, o.L1, o.L2, (lo.alpha + hi.alpha)/2; x0=lo.vector)
        x.converged || return nothing
        mid = (param=pm, alpha=x.alpha, vector=x.vector)
        imag(mid.alpha)*imag(lo.alpha) > 0 ? (lo = mid) : (hi = mid)
        abs(hi.param - lo.param) < 1e-6 && break
    end
    R, m = RM((lo.param + hi.param)/2); a = (lo.alpha + hi.alpha)/2
    zs = doppler_zeros(real(a), m/R, w/R)
    return (omega1=w, direction=String(dir), R=R, m=m, beta=m/R, alpha_real=real(a), alpha_imag=imag(a),
            wave_angle_deg=atand(m/R, real(a)), phase_speed_disk=(w/R)/hypot(real(a), m/R),
            z_crit_outer=isempty(zs) ? NaN : zs[end])
end
neutral = filter(!isnothing, pmap(b -> refine_bracket(b...), brackets))
sort!(neutral, by=x -> (x.omega1, x.R, x.m))
write_csv(joinpath(OUTDIR, "neutral_points.csv"), neutral)

open(joinpath(OUTDIR, "traveling_growth_maps_report.md"), "w") do io
    println(io, "# 行进波增长率云图与中性曲线（对照 Thomas 2024 图 2）\n")
    println(io, "频率 `omega1 = omega*R = $(OMEGA1_LIST)`（旋转系）；网格 `R=$(R_MIN):$(R_STEP):$(R_MAX)`，`m=$(M_MIN):$(M_STEP):$(M_MAX)`（仅 Re α>0，与 Thomas 2024 约定一致），共 $(length(pts)) 点；用时 $(round(time()-t0)) s。\n")
    println(io, "候选窗口 $(-IMAG_NEG) ≤ Im α ≤ $(IMAG_KEEP)。筛选：N=$(N_FIND) 全谱发现，N=$(N_MID)/$(N_REF) 局部收敛（相对差 < $(TOL)），zmax=$(ZMAX_SCAN)→$(ZMAX_CHECK) 局部检验（< $(TOL_Z)），Im α < $(IMAG_BB) 的候选做局部 Briggs–Bers，仅 plus 支、远场能量占比 < 0.1 计入。")
    println(io, @sprintf("回归：omega=0、R=285.36、beta=0.07759 通过同一筛选找回 Type-I 参考模态，|Δα|=%.2e。\n", reg_err))
    println(io, "|omega1|不稳定网格点|最大 -Im α|位置 (R, m)|中性点数|最小中性 R（m）|\n|---:|---:|---:|---|---:|---|")
    for w in OMEGA1_LIST
        u = [v for (k, v) in best if k[1] == w && v.alpha_imag < 0]
        nn = [x for x in neutral if x.omega1 == w]
        b = isempty(u) ? nothing : u[argmin([v.alpha_imag for v in u])]
        rmin = isempty(nn) ? nothing : nn[argmin([x.R for x in nn])]
        println(io, @sprintf("|%g|%d|%s|%s|%d|%s|", w, length(u), b === nothing ? "—" : @sprintf("%.4f", -b.alpha_imag),
                b === nothing ? "—" : @sprintf("(%.0f, %.0f)", b.R, b.m), length(nn),
                rmin === nothing ? "—" : @sprintf("%.1f（%.2f）", rmin.R, rmin.m)))
    end
    println(io, "\n数据：`map_modes.csv`（全部候选）、`map_least_damped.csv`（每点阻尼最小的物理模态）、`neutral_points.csv`（加密中性点，direction 表示沿 R 或沿 m 二分）。")
    println(io, "\n解释边界：云图只取每点阻尼最小的物理模态，Type-I 与 Type-II 重叠处显示较不稳定者；中性点来自相邻网格的变号区间，若两支在同一区间交替为最不稳定者，变号可能来自分支切换而非真实中性点（已用 |Δα|<0.1 排除大部分情况）。")
end
println("wrote ", OUTDIR)

#=
analyze_gap_filter_core.jl
==========================
Clean core test of the gap-filter theory (complements analyze_gap_filter.jl): inside the genuinely
uniform core 0.2 <= z <= 0.5 (Re_h = 1000, R = 330; the Bödewadt overshoot occupies z >~ 0.55) compare
  numerical  ln(|u(z_b)|/|u(z_a)|),  z_a = 0.2, z_b = 0.5  (radial velocity component, which dominates for beta << alpha)
  theory     ln(|cosh(m(zeta_b - zeta_r))| / |cosh(m(zeta_a - zeta_r))|), rigid lid at the rotor-layer edge zeta_r,
with m from (i) the inviscid formula and (ii) the exact constant-coefficient root (Re > 0 branch).
Also records the critical levels sigma(z) = beta (1 - G(z) - c) = 0 inside 0.12 < z < 0.9 (inertial regime diagnosis).
Output: results/gap_filter/core_test.csv
=#
src = read(joinpath(@__DIR__, "analyze_gap_filter.jl"), String)
eval(Meta.parseall(src[1:findfirst("Re0, R0 = 1000.0", src)[1]-1]))
b0 = cavity_baseflow(1000.0; N=300); R0 = 330.0; K0 = bary_interp(b0.z, b0.G, 0.5); zr0, zs0 = edges(b0)
s = sqrt(1000.0); za, zb = 0.2, 0.5
zq = range(0.12, 0.9, length=2001); Gq = [bary_interp(b0.z, b0.G, t) for t in zq]
rows = String[]
cases = [(a, n, c) for n in (-20, -10, 10, 20) for c in (0.70, 0.75, 0.85, 0.9, 1.0) for a in 0.02:0.04:0.62]
append!(cases, [(a, 16, NaN) for a in 0.02:0.04:0.62])          # target-like (omega_bar = 8)
for (a, n, c) in cases
    beta = n/R0; om = isnan(c) ? 8/R0 : -c*beta
    cc = isnan(c) ? -om/beta : c
    op, q = forced(b0, R0, a, n, cc)
    un = log(abs(bary_interp(op.z, comp(op, q, 1), zb))/abs(bary_interp(op.z, comp(op, q, 1), za)))
    m, sig, f = m_inviscid(a, beta, om, K0, R0)
    lr = core_roots(b0, R0, a, beta, om)[1:2]; me = lr[argmax(real.(lr))]
    th(mm) = log(abs(cosh(mm*(zb - zr0)*s))/abs(cosh(mm*(za - zr0)*s)))
    ncrit = count(k -> sign(1 - Gq[k] - cc) != sign(1 - Gq[k+1] - cc), 1:length(Gq)-1)
    push!(rows, join((a, n, cc, real(sig/f), real(m), imag(m), real(me), imag(me), un, th(m), th(me), ncrit), ","))
end
open(joinpath(@__DIR__, "results", "gap_filter", "core_test.csv"), "w") do fh
    println(fh, "alpha,n,c,sig_over_f,m_r,m_i,me_r,me_i,lnratio_num,lnratio_inviscid,lnratio_exact,n_critical_levels")
    foreach(r -> println(fh, r), rows)
end
println("wrote core_test.csv ($(length(rows)) cases)")

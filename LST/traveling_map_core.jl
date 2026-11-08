# Shared definitions for compute_traveling_growth_maps.jl (loaded on master and workers).
#
# Faster physical-mode filter for dense (R, m) maps at fixed omega1 = omega*R:
#   discovery   : one full QZ of the spatial pencil on the N_FIND grid;
#   convergence : each candidate refined locally (shift-invert, local_mode_tracking.jl)
#                 on N_MID and N_REF; accepted if |a_MID - a_FIND| < 1e-3 and
#                 |a_REF - a_MID| < TOL (relative), the reported alpha is the N_REF value;
#   zmax check  : local refinement with base-flow cut-off ZMAX_CHECK;
#   branch      : Briggs–Bers, omega -> omega + i*sigma, followed with the local solver
#                 on N_FIND, only for candidates with Im(alpha) < IMAG_BB (only these
#                 can affect the unstable region or the neutral curves).
# The full-spectrum, three-grid filter of lab_stationary_core.jl is kept unchanged there.

include(joinpath(@__DIR__, "lab_stationary_core.jl"))
include(joinpath(@__DIR__, "local_mode_tracking.jl"))

const N_FIND = envint("N_FIND", 49)   # discovery only; accuracy comes from local refinement on N_MID/N_REF
const N_MID = envint("N_MID", 89)
const N_REF = envint("N_REF", 109)
const IMAG_KEEP = envfloat("IMAG_KEEP", 0.06)   # upper bound only: damped modes beyond this are irrelevant
const IMAG_NEG = envfloat("IMAG_NEG", 0.5)      # lower bound: strongly growing modes must be kept
# Thomas (2024) counts Re(alpha) > 0 only; (m, alpha, omega) -> (-m, -conj(alpha), -omega) maps the
# Re(alpha) < 0 modes onto the same physical waves at -omega
const POSITIVE_ALPHA = get(ENV, "POSITIVE_ALPHA", "1") == "1"
const IMAG_BB = envfloat("IMAG_BB", 0.03)
const SIG_LOCAL = [1e-3*1.25^k for k in 0:28]

"""Briggs–Bers branch with the local solver; early exit once the mode is clearly in the upper half plane."""
function bb_local(R, m, omega, alpha0)
    a = alpha0; x = nothing; s_last = 0.0
    for s in SIG_LOCAL
        o = make_operator(R, m, omega; N=N_FIND, sigma=s)
        r = refine_mode(o.L0, o.L1, o.L2, a; x0=x)
        (r.converged && abs(r.alpha - a) <= 0.05) || break
        a, x, s_last = r.alpha, r.vector, s
        imag(a) > 0.2 && imag(a) > imag(alpha0) && return ("plus", imag(a), s)
    end
    dai = imag(a) - imag(alpha0)
    br = if s_last == SIG_LOCAL[end]
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

"""All zero crossings of the rotating-frame Doppler profile Re(alpha) F + beta (G-1) - omega for z < 15."""
function doppler_zeros(ar, beta, omega)
    u = ar .* FB .+ beta .* (GB .- 1.0) .- omega
    zs = Float64[]
    for j in 1:length(u)-1
        ZB[j] < 15 || break
        u[j]*u[j+1] < 0 && push!(zs, ZB[j] - u[j]*(ZB[j+1]-ZB[j])/(u[j+1]-u[j]))
    end
    zs
end

"""
Least-damped physical mode at one (R, m, omega1) point.

Candidates are examined in order of increasing Im(alpha); the loop stops at the first
physical mode (the map and the neutral curves only need that one) or when the
remaining candidates are too damped for the Briggs–Bers test (Im alpha >= IMAG_BB).
Rejected candidates examined on the way are returned as well.
"""
function map_point(omega1, R, m)
    omega = omega1/R
    op = make_operator(R, m, omega; N=N_FIND)
    modes = solve_spatial(op; shift=0.38+0im, nev=size(op.L0, 1))
    cands = [x for x in modes if isfinite(x.alpha) && abs(real(x.alpha)) <= ALPHA_MAX &&
             -IMAG_NEG <= imag(x.alpha) <= IMAG_KEEP && (!POSITIVE_ALPHA || real(x.alpha) > 0) &&
             abs(x.alpha) > 1e-4 && residual(op, x.alpha, x.vector) < RES_MAX]
    isempty(cands) && return NamedTuple[]
    opM = make_operator(R, m, omega; N=N_MID); opR = make_operator(R, m, omega; N=N_REF)
    opZ = make_operator(R, m, omega; N=N_REF, zmax=ZMAX_CHECK)
    out = NamedTuple[]; seen = ComplexF64[]
    for c in sort(cands, by=x -> imag(x.alpha))
        imag(c.alpha) >= IMAG_BB && break
        mM = refine_mode(opM.L0, opM.L1, opM.L2, c.alpha)
        (mM.converged && reldiff(mM.alpha, c.alpha) < 1e-2) || continue
        mR = refine_mode(opR.L0, opR.L1, opR.L2, mM.alpha)
        (mR.converged && reldiff(mR.alpha, mM.alpha) < TOL) || continue
        a = mR.alpha
        any(abs(a - s) < 1e-8 for s in seen) && continue
        push!(seen, a)
        mZ = refine_mode(opZ.L0, opZ.L1, opZ.L2, a)
        zok = mZ.converged && reldiff(mZ.alpha, a) < TOL_Z
        br, ai_bb, s_bb = if imag(a) < IMAG_BB
            aF = refine_mode(op.L0, op.L1, op.L2, a).alpha
            bb_local(R, m, omega, aF)
        else
            ("skipped", NaN, NaN)
        end
        zs = doppler_zeros(real(a), m/R, omega)
        zp, ff = profile_metrics(opR, mR.vector)
        push!(out, (omega1=Float64(omega1), R=Float64(R), m=Float64(m), alpha_real=real(a), alpha_imag=imag(a),
                    dN=reldiff(a, mM.alpha), dzmax=zok ? reldiff(mZ.alpha, a) : NaN, zmax_ok=zok, branch=br,
                    imag_alpha_bb=ai_bb, sigma_bb=s_bb, wave_angle_deg=atand(m/R, real(a)),
                    n_doppler_zeros=length(zs), z_crit_outer=isempty(zs) ? NaN : zs[end],
                    z_peak=zp, farfield_fraction=ff))
        map_physical(out[end]) && break
    end
    return out
end

map_physical(r) = r.zmax_ok && r.branch == "plus" && r.farfield_fraction < 0.1

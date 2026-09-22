#!/usr/bin/env julia

using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const FIT_FILE = joinpath(ROOT, "work", "results", "consistent_fold_tail_scaling",
                          "polynomial_fits.csv")
const CHAIN_FILE = joinpath(ROOT, "work", "results", "consistent_epsilon_fold_chain",
                            "recommended_epsilon_fold_table.csv")
const PROFILE_FILE = joinpath(ROOT, "work", "results", "consistent_inner_leading_matching",
                              "fixed_eta_profiles.csv")
const OUTDIR = joinpath(ROOT, "work", "results", "degenerate_inner_core_validation")
const GAMMA = 1.0

function read_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(lines[1], ',')
    rows = [Dict(header .=> split(line, ',')) for line in lines[2:end]]
    return rows
end

number(row, key) = parse(Float64, row[key])

function core_diagnostics(ro, tw; gamma=GAMMA)
    s = (2 - ro - ro^2) / 2
    omega = s + ro
    chi = 1 - gamma * (tw - 1)
    tw_theory = 1 + (1 - (omega / s)^2) / gamma
    delta_tw = tw - tw_theory
    delta_core = chi * s^2 - omega^2
    relative_core = delta_core / omega^2
    gm = (omega / sqrt(chi) - s) / ro
    return (; ro, tw, s, omega, chi, tw_theory, delta_tw,
            delta_core, relative_core, gm)
end

function polynomial_limit(epsilon, values, order)
    X = hcat((epsilon .^ k for k in 0:order)...)
    return (X \ values)[1]
end

function sign_limits(epsilon, values, half_ranges, order)
    n = length(epsilon)
    limits = Float64[]
    for mask in 0:(2^n - 1)
        perturbed = similar(values)
        for i in 1:n
            sign = ((mask >> (i - 1)) & 1) == 1 ? 1.0 : -1.0
            perturbed[i] = values[i] + sign * half_ranges[i]
        end
        push!(limits, polynomial_limit(epsilon, perturbed, order))
    end
    return limits
end

function csv_write(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            println(io, join(row, ','))
        end
    end
end

mkpath(OUTDIR)

# 1. Compare every stored extrapolation window/order with the analytic core curve.
fit_rows = read_csv(FIT_FILE)
compatibility_rows = Vector{Vector{String}}()
for window in ("all_five", "drop_epsilon_0p1"), order in (1, 2)
    ro_row = only(filter(r -> r["variable"] == "Ro" && r["window"] == window &&
                              parse(Int, r["order"]) == order, fit_rows))
    tw_row = only(filter(r -> r["variable"] == "Tw" && r["window"] == window &&
                              parse(Int, r["order"]) == order, fit_rows))
    d = core_diagnostics(number(ro_row, "limit"), number(tw_row, "limit"))
    push!(compatibility_rows, [window, string(order),
        @sprintf("%.12e", d.ro), @sprintf("%.12e", d.tw),
        @sprintf("%.12e", d.s), @sprintf("%.12e", d.omega),
        @sprintf("%.12e", d.chi), @sprintf("%.12e", d.tw_theory),
        @sprintf("%.12e", d.delta_tw), @sprintf("%.12e", d.delta_core),
        @sprintf("%.12e", d.relative_core), @sprintf("%.12e", d.gm)])
end
csv_write(joinpath(OUTDIR, "fit_window_compatibility.csv"),
    ["window", "order", "Ro_t", "Tw_t", "s_t", "omega_t", "chi_t",
     "Tw_theory", "Tw_minus_theory", "Delta_core", "Delta_core_over_omega2", "g_m"],
    compatibility_rows)

# 2. Propagate the already-recorded discretization half-ranges through quadratic fits.
chain_rows = read_csv(CHAIN_FILE)
uncertainty_rows = Vector{Vector{String}}()
for (window, keep) in (("all_five", collect(1:5)),
                       ("drop_epsilon_0p1", collect(2:5)))
    epsilon = [number(chain_rows[i], "epsilon") for i in keep]
    ro = [number(chain_rows[i], "Ro_f") for i in keep]
    tw = [number(chain_rows[i], "Tw_f") for i in keep]
    ro_half = [number(chain_rows[i], "Ro_half_range") for i in keep]
    tw_half = [number(chain_rows[i], "Tw_half_range") for i in keep]
    ro_limits = sign_limits(epsilon, ro, ro_half, 2)
    tw_limits = sign_limits(epsilon, tw, tw_half, 2)
    diagnostics = [core_diagnostics(r, t) for r in ro_limits for t in tw_limits]
    for (quantity, vals) in (
        ("Ro_t", [d.ro for d in diagnostics]),
        ("Tw_t", [d.tw for d in diagnostics]),
        ("Tw_minus_theory", [d.delta_tw for d in diagnostics]),
        ("Delta_core", [d.delta_core for d in diagnostics]),
        ("Delta_core_over_omega2", [d.relative_core for d in diagnostics]),
        ("g_m", [d.gm for d in diagnostics]))
        lo, hi = extrema(vals)
        push!(uncertainty_rows, [window, "quadratic", quantity,
            @sprintf("%.12e", lo), @sprintf("%.12e", hi),
            string(lo <= 0 <= hi), string(length(vals))])
    end
end
csv_write(joinpath(OUTDIR, "uncertainty_envelopes.csv"),
    ["window", "fit", "quantity", "minimum", "maximum", "contains_zero", "corner_count"],
    uncertainty_rows)

# The finite-epsilon points need not satisfy the limiting identity, but their defect should decay.
finite_rows = Vector{Vector{String}}()
for row in chain_rows
    epsilon = number(row, "epsilon")
    d = core_diagnostics(number(row, "Ro_f"), number(row, "Tw_f"))
    push!(finite_rows, [@sprintf("%.12e", epsilon),
        @sprintf("%.12e", d.ro), @sprintf("%.12e", d.tw),
        @sprintf("%.12e", d.tw_theory), @sprintf("%.12e", d.delta_tw),
        @sprintf("%.12e", d.delta_core), @sprintf("%.12e", d.relative_core),
        @sprintf("%.12e", d.gm)])
end
csv_write(joinpath(OUTDIR, "finite_epsilon_compatibility.csv"),
    ["epsilon", "Ro_f", "Tw_f", "Tw_theory", "Tw_minus_theory", "Delta_core",
     "Delta_core_over_omega2", "g_m"], finite_rows)

# 3. Fixed-eta profile extrapolations. This tests only the zero intercept, not first-order matching.
profile_rows = read_csv(PROFILE_FILE)
profile_fit_rows = Vector{Vector{String}}()
for eta in (1.0, 2.0, 5.0), (window, eps_cut) in (("all_five", Inf),
                                                    ("drop_epsilon_0p1", 0.05))
    selected = filter(r -> number(r, "eta") == eta && number(r, "epsilon") <= eps_cut,
                      profile_rows)
    sort!(selected, by=r -> number(r, "epsilon"), rev=true)
    epsilon = [number(r, "epsilon") for r in selected]
    field_data = (
        ("H", [number(r, "H") for r in selected]),
        ("F", [number(r, "F") for r in selected]),
        ("G", [number(r, "G") for r in selected]),
        ("Tw_minus_T", [-number(r, "T_minus_Tw") for r in selected]))
    for (field, values) in field_data
        X = hcat(ones(length(epsilon)), epsilon, epsilon.^2)
        coeff = X \ values
        residual = values - X * coeff
        last_value = values[argmin(epsilon)]
        push!(profile_fit_rows, [@sprintf("%.1f", eta), window, field,
            @sprintf("%.12e", coeff[1]), @sprintf("%.12e", coeff[2]),
            @sprintf("%.12e", coeff[3]), @sprintf("%.12e", norm(residual) / sqrt(length(values))),
            @sprintf("%.12e", last_value),
            @sprintf("%.12e", abs(coeff[1]) / max(abs(last_value), eps(Float64))),
            string(length(values))])
    end
end
csv_write(joinpath(OUTDIR, "fixed_eta_zero_intercepts.csv"),
    ["eta", "window", "field", "intercept", "linear_coefficient", "quadratic_coefficient",
     "rms_residual", "value_at_smallest_epsilon", "abs_intercept_over_last_value", "npoints"],
    profile_fit_rows)

# Last-pair effective decay exponents: a direct zero-limit test independent of a fitted intercept.
decay_rows = Vector{Vector{String}}()
for eta in (1.0, 2.0, 5.0)
    selected = filter(r -> number(r, "eta") == eta &&
                           number(r, "epsilon") in (0.01, 0.005), profile_rows)
    sort!(selected, by=r -> number(r, "epsilon"), rev=true)
    field_data = (
        ("H", [number(r, "H") for r in selected]),
        ("F", [number(r, "F") for r in selected]),
        ("G", [number(r, "G") for r in selected]),
        ("Tw_minus_T", [-number(r, "T_minus_Tw") for r in selected]))
    epsilon = [number(r, "epsilon") for r in selected]
    for (field, values) in field_data
        p_eff = log(abs(values[1] / values[2])) / log(epsilon[1] / epsilon[2])
        push!(decay_rows, [@sprintf("%.1f", eta), field,
            @sprintf("%.12e", values[1]), @sprintf("%.12e", values[2]),
            @sprintf("%.9f", p_eff)])
    end
end
csv_write(joinpath(OUTDIR, "fixed_eta_effective_decay.csv"),
    ["eta", "field", "value_at_epsilon_0p01", "value_at_epsilon_0p005", "p_effective"],
    decay_rows)

open(joinpath(OUTDIR, "scope.txt"), "w") do io
    println(io, "Algebraic/post-processing validation only; no new base-flow or continuation solve.")
    println(io, "gamma=1; input files:")
    println(io, FIT_FILE)
    println(io, CHAIN_FILE)
    println(io, PROFILE_FILE)
    println(io, "Uncertainty envelopes use all sign corners of stored numerical half-ranges;")
    println(io, "they are deterministic sensitivity envelopes, not confidence intervals.")
end

println("Core compatibility from stored extrapolations:")
for row in compatibility_rows
    @printf("  %-21s order %s: dTw=%s  Delta=%s  rel=%s  gm=%s\n",
            row[1], row[2], row[9], row[10], row[11], row[12])
end
println("\nQuadratic uncertainty envelopes:")
for row in uncertainty_rows
    if row[3] in ("Tw_minus_theory", "Delta_core", "Delta_core_over_omega2", "g_m")
        @printf("  %-21s %-24s [%s, %s], zero=%s\n",
                row[1], row[3], row[4], row[5], row[6])
    end
end
println("\nWrote results to: ", OUTDIR)

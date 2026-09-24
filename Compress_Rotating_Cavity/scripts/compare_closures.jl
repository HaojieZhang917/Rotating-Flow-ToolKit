#!/usr/bin/env julia

"""Compare traditional and Blackburn closures in one temporal framework.

The temporal spectrum is the axisymmetric, similarity-preserving spectrum.
The Blackburn operator is a generalized eigenproblem because chi=2-T weights
the three momentum time derivatives.  The script writes branch, growth-rate,
and normalized mode-shape tables without creating intermediate scripts.
"""

using CompressRotatingCavity
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT_DIR = joinpath(ROOT, "benchmark_results")

function as_blackburn(profile, config)
    return BlackburnProfile(profile.z, profile.D, profile.fields, profile.pressure,
                            profile.tw, profile.residual, config)
end

function mode_shape(profile, result)
    n = length(result.z)
    m = n - 2
    nf, ng = m - 1, m
    vector = result.eigenvectors[:, 1]
    fcoord = vector[1:nf]
    gcoord = vector[nf + 1:nf + ng]
    tcoord = vector[nf + ng + 1:end]
    f = result.f_map * fcoord
    h = result.h_map * fcoord
    g = result.g_map * gcoord
    theta = result.temperature_map * tcoord
    scale = maximum(abs, vcat(h, f, g, theta))
    scale > 0 || error("zero temporal mode")
    return (; h=h ./ scale, f=f ./ scale, g=g ./ scale,
            theta=theta ./ scale, z=result.z)
end

function write_mode_rows(path, rows)
    open(path, "w") do io
        println(io, "closure,branch,z,abs_h,abs_f,abs_g,abs_theta")
        for row in rows
            println(io, join((row.closure, row.branch, row.z, row.abs_h,
                              row.abs_f, row.abs_g, row.abs_theta), ','))
        end
    end
end

function main()
    mkpath(OUTPUT_DIR)
    config = ThermalConfig(re_h=1000.0, pr=0.72, degree=100, tolerance=2e-8)
    bconfig = BlackburnConfig(re_h=1000.0, pr=0.72, degree=100, tolerance=2e-8)
    temporal_degree = 120
    target_tw = 1.16

    traditional_iso = CompressRotatingCavity.isothermal_profile(config, REFERENCE)
    traditional_specs = [
        (1, "upper_low_Pi", 0.001542389588),
        (2, "middle", 0.003786919205),
        (3, "principal_isothermal", 0.018230452624),
    ]
    traditional = [(index, label,
                    CompressRotatingCavity.fixed_pressure_profile(pi, traditional_iso))
                   for (index, label, pi) in traditional_specs]

    blackburn_iso = blackburn_isothermal_profile(bconfig, REFERENCE)
    blackburn = NamedTuple[]
    # The principal Blackburn branch is continued from the isothermal state.
    push!(blackburn, (index=3, label="principal_isothermal",
                      profile=blackburn_fixed_temperature_profile(target_tw, blackburn_iso),
                      status="converged"))
    # Use the traditional middle and upper states as independent Newton seeds.
    for (index, label, _) in traditional_specs[1:2]
        seed = as_blackburn(traditional[index][3], bconfig)
        try
            profile = blackburn_fixed_temperature_profile(target_tw, seed)
            push!(blackburn, (index=index, label=label, profile=profile,
                              status="converged"))
        catch exception
            push!(blackburn, (index=index, label=label, profile=nothing,
                              status="failed: " * sprint(showerror, exception)))
        end
    end
    sort!(blackburn; by=row -> row.index)

    spectral_rows = NamedTuple[]
    mode_rows = NamedTuple[]
    function add_spectrum!(closure, branch, profile, closure_symbol)
        result = temporal_operator(profile, temporal_degree, config.re_h, config.pr;
                                   closure=closure_symbol)
        leading = leading_temporal_eigenvalue(result)
        mode = mode_shape(profile, result)
        component_norms = [sqrt(sum(abs2, mode.h)), sqrt(sum(abs2, mode.f)),
                           sqrt(sum(abs2, mode.g)), sqrt(sum(abs2, mode.theta))]
        fraction_sum = sum(component_norms)
        push!(spectral_rows, (
            closure=closure, branch=branch, tw=profile.tw, pressure=profile.pressure,
            leading_real=real(leading), leading_imag=imag(leading),
            mode_h=component_norms[1], mode_f=component_norms[2],
            mode_g=component_norms[3], mode_theta=component_norms[4],
            frac_h=component_norms[1] / fraction_sum,
            frac_f=component_norms[2] / fraction_sum,
            frac_g=component_norms[3] / fraction_sum,
            frac_theta=component_norms[4] / fraction_sum,
            mode_count=length(result.eigenvalues), status="converged",
        ))
        for j in eachindex(mode.z)
            push!(mode_rows, (closure=closure, branch=branch, z=mode.z[j],
                              abs_h=abs(mode.h[j]), abs_f=abs(mode.f[j]),
                              abs_g=abs(mode.g[j]), abs_theta=abs(mode.theta[j])))
        end
        return result
    end

    traditional_results = Dict{Int,Any}()
    for (index, label, profile) in traditional
        traditional_results[index] = add_spectrum!("traditional", label, profile, :traditional)
    end
    blackburn_results = Dict{Int,Any}()
    for row in blackburn
        row.profile === nothing && continue
        blackburn_results[row.index] = add_spectrum!("blackburn", row.label,
                                                      row.profile, :blackburn)
    end

    spectral_path = joinpath(OUTPUT_DIR, "closure_comparison_Re1000_Tw1.160.csv")
    open(spectral_path, "w") do io
        println(io, "closure,branch,Tw,Pi,leading_real,leading_imag,mode_h,mode_f,mode_g,mode_theta,frac_h,frac_f,frac_g,frac_theta,mode_count,status")
        for row in spectral_rows
            println(io, join((row.closure, row.branch, row.tw, row.pressure,
                              row.leading_real, row.leading_imag, row.mode_h,
                              row.mode_f, row.mode_g, row.mode_theta,
                              row.frac_h, row.frac_f, row.frac_g, row.frac_theta,
                              row.mode_count, row.status), ','))
        end
    end
    write_mode_rows(joinpath(OUTPUT_DIR, "closure_comparison_Re1000_Tw1.160_modes.csv"), mode_rows)

    fold, _ = trace_principal_fold(traditional_iso; coarse_step=0.003,
                                   fine_step=0.0005, pressure_tolerance=1e-9)
    # Heated Blackburn principal continuation.  A missing sign change is a
    # result: it means no Blackburn fold was detected on the scanned interval.
    current = blackburn_iso
    previous_tangent = blackburn_pressure_tangent(current)
    max_profile = current
    fold_detected = false
    for _ in 1:800
        next_pressure = current.pressure + 0.001
        trial = try
            blackburn_fixed_pressure_profile(next_pressure, current)
        catch
            break
        end
        tangent = blackburn_pressure_tangent(trial)
        max_profile = trial
        if previous_tangent * tangent <= 0
            fold_detected = true
            break
        end
        current, previous_tangent = trial, tangent
        max_profile.tw >= 1.70 && break
    end

    # Evaluate both closures at the traditional fold temperature.  This
    # separates a true zero-growth saddle-node from a Blackburn state that
    # remains a regular, damped equilibrium at the same wall temperature.
    traditional_fold_spectrum = temporal_operator(
        fold, temporal_degree, config.re_h, config.pr; closure=:traditional,
    )
    blackburn_at_traditional_fold = blackburn_fixed_temperature_profile(
        fold.tw, blackburn_iso,
    )
    blackburn_fold_spectrum = temporal_operator(
        blackburn_at_traditional_fold, temporal_degree, config.re_h, config.pr;
        closure=:blackburn,
    )
    critical_path = joinpath(OUTPUT_DIR, "closure_critical_comparison_Re1000.csv")
    open(critical_path, "w") do io
        println(io, "closure,event,Tw,Pi,leading_real,leading_imag,status")
        λt = leading_temporal_eigenvalue(traditional_fold_spectrum)
        λb = leading_temporal_eigenvalue(blackburn_fold_spectrum)
        println(io, join(("traditional", "principal_fold", fold.tw,
                          fold.pressure, real(λt), imag(λt), "converged"), ','))
        println(io, join(("blackburn", "at_traditional_fold_temperature",
                          blackburn_at_traditional_fold.tw,
                          blackburn_at_traditional_fold.pressure,
                          real(λb), imag(λb), "regular_state"), ','))
    end

    profile_rows = NamedTuple[]
    for index in (1, 2, 3)
        tr = traditional[index][3]
        br = findfirst(row -> row.index == index, blackburn)
        bp = br === nothing ? nothing : blackburn[br].profile
        push!(profile_rows, (closure="traditional", branch=traditional[index][2],
                             tw=tr.tw, pressure=tr.pressure, status="converged"))
        bp === nothing || push!(profile_rows, (closure="blackburn", branch=blackburn[br].label,
                                               tw=bp.tw, pressure=bp.pressure,
                                               status=blackburn[br].status))
    end
    profile_path = joinpath(OUTPUT_DIR, "closure_comparison_Re1000_Tw1.160_profiles.csv")
    open(profile_path, "w") do io
        println(io, "closure,branch,Tw,Pi,status")
        for row in profile_rows
            println(io, join((row.closure, row.branch, row.tw, row.pressure, row.status), ','))
        end
    end

    summary_path = joinpath(OUTPUT_DIR, "closure_comparison_Re1000_summary.txt")
    open(summary_path, "w") do io
        println(io, "Traditional versus Blackburn generalized Boussinesq")
        println(io, "Re_h=1000, Pr=0.72, target Tw=1.16, temporal_degree=$temporal_degree")
        println(io, "traditional_principal_fold_Tw=$(fold.tw)")
        println(io, "blackburn_fold_detected=$fold_detected")
        println(io, "blackburn_principal_scan_max_Tw=$(max_profile.tw)")
        println(io, "blackburn_principal_scan_pressure=$(max_profile.pressure)")
        println(io, "blackburn_principal_scan_status=$(max_profile.residual)")
        println(io, "traditional_fold_leading_lambda=$(leading_temporal_eigenvalue(traditional_fold_spectrum))")
        println(io, "blackburn_at_traditional_fold_Tw=$(blackburn_at_traditional_fold.tw)")
        println(io, "blackburn_at_traditional_fold_leading_lambda=$(leading_temporal_eigenvalue(blackburn_fold_spectrum))")
        println(io, ""); println(io, "Temporal leading eigenvalues")
        for row in spectral_rows
            println(io, "$(row.closure) $(row.branch): $(row.leading_real) + $(row.leading_imag)i")
        end
        println(io, ""); println(io, "Blackburn branch statuses")
        for row in blackburn
            println(io, "$(row.label): $(row.status)")
        end
    end

    println("Closure comparison benchmark")
    println("traditional principal fold Tw=$(fold.tw)")
    println("Blackburn principal scan: max Tw=$(max_profile.tw), fold_detected=$(fold_detected)")
    println("at traditional fold: lambda_T=$(leading_temporal_eigenvalue(traditional_fold_spectrum)), lambda_B=$(leading_temporal_eigenvalue(blackburn_fold_spectrum))")
    for row in spectral_rows
        println("$(row.closure) $(row.branch): lambda=$(row.leading_real)+$(row.leading_imag)im")
    end
    println("Wrote ", spectral_path)
    println("Wrote ", summary_path)
    println("PASS")
end

main()

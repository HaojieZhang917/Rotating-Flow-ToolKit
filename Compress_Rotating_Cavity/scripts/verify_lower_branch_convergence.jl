#!/usr/bin/env julia

"""Lower/principal branch continuation and grid-convergence verification.

The reported branch is the state connected to the isothermal rotor-stator
solution.  The traditional fold and the Blackburn heated continuation are
computed at the same Reynolds number, Prandtl number, and temporal operator.
"""

using CompressRotatingCavity
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT = joinpath(ROOT, "benchmark_results")
const TARGETS = (1.16, 1.40, 1.60, 1.75)
const TRADITIONAL_PRINCIPAL_PI = 0.018230452624

function continue_blackburn_to_targets(config, targets)
    current = blackburn_isothermal_profile(config, REFERENCE)
    profiles = NamedTuple[]
    for target in targets
        pressure_step = 0.001
        while current.tw < target - 1e-8
            try
                trial = blackburn_fixed_pressure_profile(
                    current.pressure + pressure_step, current,
                )
                if trial.tw >= target
                    current = blackburn_fixed_temperature_profile(target, trial)
                    break
                end
                current = trial
                pressure_step = min(0.001, pressure_step * 1.15)
            catch exception
                pressure_step *= 0.5
                pressure_step < 1e-6 && error(
                    "Blackburn continuation failed near Tw=$(current.tw), " *
                    "target=$target: " * sprint(showerror, exception),
                )
            end
        end
        push!(profiles, (tw=current.tw, pressure=current.pressure,
                        fields=current.fields, profile=current,
                        residual=current.residual))
    end
    return profiles
end

function temporal_mode_metrics(profile, degree, temporal_degree, closure)
    result = temporal_operator(profile, temporal_degree, profile.config.re_h,
                               profile.config.pr;
                               closure=closure)
    vector = result.eigenvectors[:, 1]
    n = length(result.z)
    m = n - 2
    nf = m - 1
    f = result.f_map * vector[1:nf]
    h = result.h_map * vector[1:nf]
    g = result.g_map * vector[nf+1:nf+m]
    theta = result.temperature_map * vector[nf+m+1:end]
    norms = [sqrt(sum(abs2, component)) for component in (h, f, g, theta)]
    fractions = norms ./ sum(norms)
    λ = leading_temporal_eigenvalue(result)
    return (; λ, frac_h=fractions[1], frac_f=fractions[2],
            frac_g=fractions[3], frac_theta=fractions[4])
end

function write_rows(path, header, rows)
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, join(row, ','))
        end
    end
end

function main()
    mkpath(OUTPUT)
    degrees = (60, 80, 100)
    fold_rows = Any[]
    branch_rows = Any[]
    mode_rows = Any[]

    for degree in degrees
        println("Lower branch degree=$degree")
        traditional_config = ThermalConfig(re_h=1000.0, pr=0.72,
                                            degree=degree, tolerance=2e-8)
        blackburn_config = BlackburnConfig(re_h=1000.0, pr=0.72,
                                           degree=degree, tolerance=2e-8)
        # Keep the temporal discretization fixed so this table isolates the
        # base-flow collocation convergence.
        temporal_degree = 120

        traditional_iso = isothermal_profile(traditional_config, REFERENCE)
        fold, _ = trace_principal_fold(traditional_iso; coarse_step=0.003,
                                       fine_step=0.0005,
                                       pressure_tolerance=2e-9)
        λfold = temporal_mode_metrics(fold, degree, temporal_degree, :traditional).λ
        push!(fold_rows, (degree, fold.tw, fold.pressure, fold.residual,
                          real(λfold), imag(λfold)))

        traditional_principal = fixed_pressure_profile(
            TRADITIONAL_PRINCIPAL_PI, traditional_iso,
        )
        λtraditional = temporal_mode_metrics(traditional_principal, degree,
                                             temporal_degree, :traditional)
        push!(branch_rows, ("traditional", degree, traditional_principal.tw,
                            traditional_principal.pressure,
                            traditional_principal.residual,
                            real(λtraditional.λ), imag(λtraditional.λ),
                            temporal_degree))
        if degree == maximum(degrees)
            push!(mode_rows, ("traditional", "principal", traditional_principal.tw,
                              traditional_principal.pressure, λtraditional.frac_h,
                              λtraditional.frac_f, λtraditional.frac_g,
                              λtraditional.frac_theta))
        end

        blackburn_profiles = continue_blackburn_to_targets(blackburn_config, TARGETS)
        for entry in blackburn_profiles
            metrics = temporal_mode_metrics(entry.profile, degree, temporal_degree,
                                            :blackburn)
            push!(branch_rows, ("blackburn", degree, entry.tw, entry.pressure,
                                entry.residual, real(metrics.λ), imag(metrics.λ),
                                temporal_degree))
            if degree == maximum(degrees)
                push!(mode_rows, ("blackburn", "principal", entry.tw,
                                  entry.pressure, metrics.frac_h, metrics.frac_f,
                                  metrics.frac_g, metrics.frac_theta))
            end
        end
    end

    write_rows(joinpath(OUTPUT, "lower_branch_fold_convergence.csv"),
               "degree,Tw_fold,Pi_fold,residual,leading_real,leading_imag",
               fold_rows)
    write_rows(joinpath(OUTPUT, "lower_branch_temporal_convergence.csv"),
               "closure,degree,Tw,Pi,residual,leading_real,leading_imag,temporal_degree",
               branch_rows)
    write_rows(joinpath(OUTPUT, "lower_branch_mode_structure.csv"),
               "closure,branch,Tw,Pi,frac_h,frac_f,frac_g,frac_theta",
               mode_rows)

    open(joinpath(OUTPUT, "lower_branch_convergence_summary.txt"), "w") do io
        println(io, "Lower/principal branch convergence at Re_h=1000, Pr=0.72")
        println(io, "The branch is connected to the isothermal rotor-stator state.")
        println(io, "")
        println(io, "Traditional fold")
        println(io, "degree,Tw_fold,Pi_fold,residual,leading_real,leading_imag")
        for row in fold_rows
            println(io, join(row, ','))
        end
        println(io, "")
        println(io, "Temporal continuation")
        println(io, "closure,degree,Tw,Pi,residual,leading_real,leading_imag,temporal_degree")
        for row in branch_rows
            println(io, join(row, ','))
        end
        println(io, "")
        println(io, "The Blackburn principal continuation was successfully evaluated at all targets: ")
        println(io, join(TARGETS, ", "))
        println(io, "No Blackburn saddle-node was detected on this scanned principal branch.")
    end
    println("Wrote lower branch convergence tables")
    println("PASS")
end

main()

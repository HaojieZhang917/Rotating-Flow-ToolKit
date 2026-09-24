#!/usr/bin/env julia

"""Scan the traditional Boussinesq temporal LST spectrum along the lower branch.

The base flow is solved first at each prescribed wall temperature.  The local
perturbation problem is then assembled at fixed `(Re_h,R,alpha,beta)`, so only
the basic state changes during the scan.  The two anchors are the isothermal
Type-I and Type-II points reported in Table IV of the manuscript for
`Re_h=1000` and `a_s=0`.
"""

using LinearAlgebra
using NonlinearEigenproblems
using Printf
using Test

const ROOT = normpath(joinpath(@__DIR__, ".."))
using CompressRotatingCavity

const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT_DIR = joinpath(ROOT, "benchmark_results")

struct ModeAnchor
    label::String
    radius::Float64
    alpha::Float64
    beta::Float64
    omega_reference::Float64
end

const ANCHORS = ModeAnchor[
    # Manuscript Table IV, a_s=0, Type I.
    ModeAnchor("Type_I", 283.21, 0.403, 0.0824, 0.0057),
    # Manuscript Table IV, a_s=0, Type II.
    ModeAnchor("Type_II", 91.08, 0.256, -0.1135, 0.1013),
]

const TEMPERATURES = vcat(collect(1.00:0.01:1.16), [1.165, 1.166, 1.167])
const MAX_PHYSICAL_FREQUENCY = 1.0e6

function finite_spectrum(cof, anchor)
    A, B = temporal_matrices(cof; alpha=anchor.alpha, beta=anchor.beta)
    values, vectors = eigen(A, B)
    # The pressure row makes B singular. LAPACK represents the associated
    # infinite eigenvalues as finite numbers of order 1e17, so remove that
    # numerical infinity before ordering the physical spectrum.
    finite = findall(i -> isfinite(values[i]) &&
                            abs(values[i]) < MAX_PHYSICAL_FREQUENCY,
                     eachindex(values))
    isempty(finite) && error("no finite temporal eigenvalues for $(anchor.label)")
    # With exp(-i*omega*t), temporal growth is Im(omega).  Keep the most
    # amplified finite mode first; the paper's branch is tracked separately.
    order = sort(finite; by=i -> (-imag(values[i]), abs(real(values[i]))))
    return values[order], vectors[:, order]
end

function component_fractions(cof, values, vectors, index)
    mode = expand_mode(cof, vectors[:, index])
    norms = [norm(mode.u), norm(mode.v), norm(mode.w), norm(mode.theta)]
    total = sum(norms)
    total > 0 || error("zero perturbation mode")
    return norms ./ total
end

function tracked_index(values, target)
    argmin(abs.(values .- target))
end

function spatial_alpha(cof, anchor, previous)
    L0, L1, L2 = spatial_matrices(cof; beta=anchor.beta,
                                   omega=anchor.omega_reference)
    values, _ = iar(PEP([L0, L1, L2]); σ=previous, neigs=1,
                    maxit=500, tol=1e-10)
    isempty(values) && error("spatial PEP returned no eigenvalue for $(anchor.label)")
    return ComplexF64(values[1])
end

function continue_temperature(seed, target; initial_step=0.002, minimum_step=1e-5)
    current = seed
    step = initial_step
    while current.tw < target - 1e-11
        next_tw = min(target, current.tw + step)
        try
            current = fixed_temperature_profile(next_tw, current)
            step = min(initial_step, 1.5 * step)
        catch exception
            step *= 0.5
            step >= minimum_step || rethrow(exception)
        end
    end
    return current
end

function main()
    mkpath(OUTPUT_DIR)
    config = ThermalConfig(re_h=1000.0, pr=0.72, degree=100, tolerance=2e-8)
    current = isothermal_profile(config, REFERENCE)

    spectra_rows = NamedTuple[]
    spatial_rows = NamedTuple[]
    tracked = Dict(anchor.label => ComplexF64(anchor.omega_reference) for anchor in ANCHORS)
    spatial_tracked = Dict(anchor.label => ComplexF64(anchor.alpha) for anchor in ANCHORS)

    for tw in TEMPERATURES
        profile = tw == 1.0 ? current : continue_temperature(current, tw)
        current = profile
        @printf("Tw=%.6f, Pi=%+.10e, residual=%.3e\n",
                profile.tw, profile.pressure, profile.residual)
        @test abs(profile.tw - tw) < 2e-8
        @test profile.residual < 2e-7

        for anchor in ANCHORS
            cof = lst_coefficients(profile; radius=anchor.radius)
            values, vectors = finite_spectrum(cof, anchor)
            previous = tracked[anchor.label]
            index = tracked_index(values, previous)
            tracked[anchor.label] = values[index]
            fractions = component_fractions(cof, values, vectors, index)
            most_amplified = values[1]
            push!(spectra_rows, (
                closure="traditional",
                mode=anchor.label,
                Tw=profile.tw,
                pressure=profile.pressure,
                radius=anchor.radius,
                alpha=anchor.alpha,
                beta=anchor.beta,
                max_growth_real=real(most_amplified),
                max_growth_imag=imag(most_amplified),
                tracked_real=real(values[index]),
                tracked_imag=imag(values[index]),
                tracked_distance=abs(values[index] - previous),
                u_fraction=fractions[1],
                v_fraction=fractions[2],
                w_fraction=fractions[3],
                theta_fraction=fractions[4],
                mode_rank=index,
            ))

            previous_alpha = spatial_tracked[anchor.label]
            alpha = spatial_alpha(cof, anchor, previous_alpha)
            spatial_tracked[anchor.label] = alpha
            push!(spatial_rows, (
                closure="traditional",
                mode=anchor.label,
                Tw=profile.tw,
                pressure=profile.pressure,
                radius=anchor.radius,
                omega=anchor.omega_reference,
                beta=anchor.beta,
                alpha_real=real(alpha),
                alpha_imag=imag(alpha),
                alpha_distance=abs(alpha - previous_alpha),
            ))
        end
    end

    spatial_csv_path = joinpath(OUTPUT_DIR,
                                "traditional_lst_spatial_temperature_scan_Re1000.csv")
    open(spatial_csv_path, "w") do io
        println(io, "closure,mode,Tw,Pi,R,omega,beta,alpha_real,alpha_imag,alpha_distance")
        for row in spatial_rows
            println(io, join((row.closure, row.mode, row.Tw, row.pressure,
                              row.radius, row.omega, row.beta,
                              row.alpha_real, row.alpha_imag,
                              row.alpha_distance), ','))
        end
    end

    csv_path = joinpath(OUTPUT_DIR, "traditional_lst_temperature_scan_Re1000.csv")
    open(csv_path, "w") do io
        println(io, "closure,mode,Tw,Pi,R,alpha,beta,max_growth_real,max_growth_imag,tracked_real,tracked_imag,tracked_distance,u_fraction,v_fraction,w_fraction,theta_fraction,mode_rank")
        for row in spectra_rows
            println(io, join((row.closure, row.mode, row.Tw, row.pressure,
                              row.radius, row.alpha, row.beta,
                              row.max_growth_real, row.max_growth_imag,
                              row.tracked_real, row.tracked_imag,
                              row.tracked_distance, row.u_fraction,
                              row.v_fraction, row.w_fraction,
                              row.theta_fraction, row.mode_rank), ','))
        end
    end

    summary_path = joinpath(OUTPUT_DIR, "traditional_lst_temperature_scan_Re1000.txt")
    open(summary_path, "w") do io
        println(io, "Traditional Boussinesq local temporal LST temperature scan")
        println(io, "Re_h=1000, Pr=0.72, base-flow degree=100")
        println(io, "lower/principal branch, Tw=1.0 to 1.167")
        println(io, "fixed perturbation parameters are those in manuscript Table IV, a_s=0")
        for anchor in ANCHORS
            rows = filter(row -> row.mode == anchor.label, spectra_rows)
            first_row, last_row = first(rows), last(rows)
            println(io)
            println(io, "$(anchor.label): R=$(anchor.radius), alpha=$(anchor.alpha), beta=$(anchor.beta)")
            @printf(io, "  Tw=1.000: tracked_omega=%+.12e%+.12ei, leading=%+.12e%+.12ei\n",
                    first_row.tracked_real, first_row.tracked_imag,
                    first_row.max_growth_real, first_row.max_growth_imag)
            @printf(io, "  Tw=%.3f: tracked_omega=%+.12e%+.12ei, leading=%+.12e%+.12ei\n",
                    last_row.Tw, last_row.tracked_real, last_row.tracked_imag,
                    last_row.max_growth_real, last_row.max_growth_imag)
            @printf(io, "  max tracked growth over scan = %.12e\n",
                    maximum(row.tracked_imag for row in rows))
            @printf(io, "  max physical growth over scan = %.12e\n",
                    maximum(row.max_growth_imag for row in rows))
        end
        println(io)
        println(io, "Spatial alpha scan at fixed omega and beta")
        for anchor in ANCHORS
            rows = filter(row -> row.mode == anchor.label, spatial_rows)
            first_row, last_row = first(rows), last(rows)
            println(io, "$(anchor.label): R=$(anchor.radius), omega=$(anchor.omega_reference), beta=$(anchor.beta)")
            @printf(io, "  Tw=1.000: alpha=%+.12e%+.12ei\n",
                    first_row.alpha_real, first_row.alpha_imag)
            @printf(io, "  Tw=%.3f: alpha=%+.12e%+.12ei\n",
                    last_row.Tw, last_row.alpha_real, last_row.alpha_imag)
            @printf(io, "  delta_alpha=%+.12e%+.12ei\n",
                    last_row.alpha_real - first_row.alpha_real,
                    last_row.alpha_imag - first_row.alpha_imag)
        end
    end

    println("Wrote ", csv_path)
    println("Wrote ", spatial_csv_path)
    println("Wrote ", summary_path)
    println("PASS")
end

main()

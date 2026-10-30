#!/usr/bin/env julia

"""Independent numerical checks for the non-similar single-disk radial march.

Run with `julia --project=. scripts/verify_single_disk_radial.jl`.
Only Julia standard libraries and the directly included radial module are used.
Results are written separately from the main attraction benchmark.
"""

using LinearAlgebra
using Test

include(joinpath(@__DIR__, "..", "src", "SingleDiskRadial.jl"))
using .SingleDiskRadial

const OUTPUT = normpath(joinpath(@__DIR__, "..", "benchmark_results", "single_disk_radial"))
const TARGET_S = log(4.0)

function check!(rows, label, value, tolerance, passed; details="")
    push!(rows, (; check=label, value, tolerance, pass=Bool(passed), details))
    println(label, ": ", value, " (", passed ? "PASS" : "FAIL", ")")
    return passed
end

function check_small!(rows, label, value, tolerance; details="")
    return check!(rows, label, value, tolerance,
                  isfinite(value) && abs(value) <= tolerance; details)
end

function csv_cell(value)
    text = string(value)
    return occursin(r"[,\"\n\r]", text) ? "\"" * replace(text, "\"" => "\"\"") * "\"" : text
end

function save_results(rows)
    mkpath(OUTPUT)
    open(joinpath(OUTPUT, "verification.csv"), "w") do io
        println(io, "check,value,tolerance,pass,details")
        for row in rows
            println(io, join(csv_cell.((row.check, row.value, row.tolerance,
                                      row.pass, row.details)), ','))
        end
    end
    open(joinpath(OUTPUT, "verification.txt"), "w") do io
        println(io, "Single-disk non-similar steady boundary-layer verification")
        println(io, "Coordinate s=log(r/r0); endpoint r/r0=4 unless stated.")
        println(io, "These checks validate this restricted steady axisymmetric model, not full NS stability.")
        println(io, "H at the inlet is uncomputed; H checks use downstream stations only.")
        println(io, "Checks: ", length(rows), "; passed: ", count(row -> row.pass, rows))
        for row in rows
            println(io, row.pass ? "PASS " : "FAIL ", row.check, " = ", row.value,
                    "; criterion/tolerance: ", row.tolerance)
        end
    end
end

function run_case(; degree=80, eta_max=20.0, ds=0.02, case=:large_radial)
    config = RadialConfig(; degree, eta_max, ds, s_end=TARGET_S)
    reference = reference_profile(config)
    result = march(reference, inlet_profile(reference; case))
    result.status == :completed || error(
        "verification case $case did not complete: $(result.status): $(result.message)")
    return result
end

"""Compare actual final physical-eta profiles on a shared sampling interval."""
function sampled_final(result, points; include_h=false)
    F = interpolate_profile(result.reference.eta, result.F[:, end], points)
    G = interpolate_profile(result.reference.eta, result.G[:, end], points)
    if include_h
        H = interpolate_profile(result.reference.eta, result.H[:, end], points)
        return vcat(F, G, H)
    end
    return vcat(F, G)
end

relative_difference(a, b) = norm(a .- b) / max(norm(b), eps(Float64))

function main()
    BLAS.set_num_threads(1)
    rows = NamedTuple[]
    try
        config = RadialConfig(degree=80, eta_max=20.0, ds=0.02, s_end=TARGET_S)
        reference = reference_profile(config)
        Fp0, Gp0 = (reference.D * reference.F)[1], (reference.D * reference.G)[1]
        # Independently obtained solve_bvp values for the same finite eta domain.
        check_small!(rows, "reference_Fp0_eta20_absolute_error", Fp0-0.51023261617, 2e-8)
        check_small!(rows, "reference_Gp0_eta20_absolute_error", Gp0+0.61592201262, 2e-8)
        check_small!(rows, "reference_Htop_eta20_absolute_error", reference.H[end]+0.8844734148, 2e-8)

        # Check the genuinely nonsteady Jacobian, with nonzero a and nonzero histories.
        eta = reference.eta
        n = length(eta)
        direction = vcat(cos.(0.6 .* eta), sin.(0.8 .* eta) .+ 0.1,
                         cos.(0.4 .* eta) .- 0.2)
        direction ./= norm(direction)
        state = reference.state .+ 0.015 .* direction
        a = 17.5
        history_f = -a .* reference.F .+ 0.07 .* sin.(pi .* eta ./ config.eta_max)
        history_g = -a .* reference.G .- 0.11 .* sin.(2pi .* eta ./ config.eta_max)
        _, J = radial_system(state, reference; a, history_f, history_g)
        step = 1e-4
        plus, _ = radial_system(state .+ step .* direction, reference;
                                a, history_f, history_g, jacobian=false)
        minus, _ = radial_system(state .- step .* direction, reference;
                                 a, history_f, history_g, jacobian=false)
        fd = (plus .- minus) ./ (2step)
        jacobian_error = norm(fd .- J * direction) / max(norm(J * direction), 1.0)
        check_small!(rows, "nonzero_history_Jacobian_directional_relative_error",
                     jacobian_error, 2e-7;
                     details="centered finite difference; nonzero inlet boundary directions included")

        exact = march(reference, inlet_profile(reference; case=:similarity))
        check!(rows, "similarity_inlet_status", exact.status, "completed", exact.status == :completed)
        check_small!(rows, "similarity_inlet_max_FG_drift",
                     maximum(metric.error_fg for metric in exact.metrics), 2e-8)
        check!(rows, "inlet_H_is_explicitly_uncomputed", all(isnan, exact.H[:, 1]),
               "true", all(isnan, exact.H[:, 1]))

        large = run_case()
        check!(rows, "large_radial_status", large.status, "completed", large.status == :completed)
        finite_fields = all(isfinite, large.F) && all(isfinite, large.G) &&
                        all(isfinite, large.H[:, 2:end])
        check!(rows, "large_radial_downstream_fields_finite", finite_fields, "true", finite_fields)
        minimum_f = minimum(metric.min_F for metric in large.metrics)
        check!(rows, "large_radial_minimum_F", minimum_f,
               ">= -$(config.reverse_tolerance)", minimum_f >= -config.reverse_tolerance)
        contraction = large.metrics[end].error_fg / large.metrics[1].error_fg
        check_small!(rows, "large_radial_FG_error_contraction_at_ratio4", contraction, 0.2;
                     details="sampled attraction for this inlet; not a general theorem")
        check_small!(rows, "large_radial_integral_continuity_max_abs",
                     maximum(abs(metric.continuity_balance) for metric in large.metrics[2:end]),
                     2e-7)

        negative_inlet = inlet_profile(reference; case=:similarity)
        bad_f = copy(negative_inlet.F)
        bad_f[argmax(reference.F)] = -0.05
        reversed = march(reference, merge(negative_inlet, (; F=bad_f)))
        check!(rows, "negative_F_inlet_reports_reverse_flow", reversed.status,
               "reverse_flow", reversed.status == :reverse_flow)
        check!(rows, "negative_F_inlet_stops_before_marching", length(reversed.s),
               "1 station", length(reversed.s) == 1)

        # Radial convergence uses the same degree/domain/inlet and the same endpoint.
        radial = [run_case(degree=64, ds=ds) for ds in (0.04, 0.02, 0.01)]
        radial_profiles = [vcat(result.F[:, end], result.G[:, end]) for result in radial]
        d_coarse = norm(radial_profiles[1] .- radial_profiles[2])
        d_fine = norm(radial_profiles[2] .- radial_profiles[3])
        order = log2(d_coarse / d_fine)
        step_ratios = [radial[1].ds / radial[2].ds, radial[2].ds / radial[3].ds]
        check!(rows, "BDF2_radial_effective_steps", join((result.ds for result in radial), ';'),
               "positive steps; successive step ratios in [1.9,2.1]",
               all(result.ds > 0 for result in radial) && all(1.9 .<= step_ratios .<= 2.1);
               details="steps are adjusted slightly to hit exactly r/r0=4")
        check!(rows, "BDF2_radial_observed_log2_order", order, "1.6 <= p <= 2.5",
               isfinite(order) && 1.6 <= order <= 2.5;
               details="final F,G differences; BE is used only at the first station")
        check!(rows, "BDF2_radial_finer_step_reduces_difference", d_fine / d_coarse,
               "< 1", isfinite(d_fine/d_coarse) && d_fine < d_coarse)

        # Spectral-grid and outer-domain checks are separate from ds convergence.
        points_grid = collect(range(0.0, 20.0; length=401))
        grid_results = [run_case(degree=N, case=:small_shape) for N in (48, 64, 80)]
        grid_profiles = [sampled_final(result, points_grid; include_h=true) for result in grid_results]
        dg_coarse = relative_difference(grid_profiles[1], grid_profiles[2])
        dg_fine = relative_difference(grid_profiles[2], grid_profiles[3])
        check_small!(rows, "eta_grid_N48_N64_relative_profile_difference", dg_coarse, 2e-5)
        check_small!(rows, "eta_grid_N64_N80_relative_profile_difference", dg_fine, 2e-7)
        check!(rows, "eta_grid_refinement_decreases_difference", dg_fine,
               "<= max(0.8*coarse_difference;1e-8)",
               dg_fine <= max(0.8dg_coarse, 1e-8);
               details="absolute floor allows nonlinear solver roundoff")

        points_domain = collect(range(0.0, 12.0; length=301))
        domain_results = [run_case(degree=80, eta_max=L, case=:small_shape) for L in (16.0, 20.0, 24.0)]
        domain_profiles = [sampled_final(result, points_domain; include_h=true) for result in domain_results]
        dl_coarse = relative_difference(domain_profiles[1], domain_profiles[2])
        dl_fine = relative_difference(domain_profiles[2], domain_profiles[3])
        check_small!(rows, "domain_eta16_eta20_inner_profile_difference", dl_coarse, 1e-4;
                     details="F,G,H compared on the shared physical interval 0<=eta<=12")
        check_small!(rows, "domain_eta20_eta24_inner_profile_difference", dl_fine, 3e-6)
        check!(rows, "domain_extension_decreases_difference", dl_fine,
               "< eta16_eta20 difference", dl_fine < dl_coarse)
        check_small!(rows, "domain_eta20_eta24_reference_Htop_difference",
                     domain_results[2].reference.H[end]-domain_results[3].reference.H[end], 2e-6;
                     details="entrainment is more sensitive to eta truncation than wall gradients")
    finally
        save_results(rows)
    end

    @testset "Single-disk radial verification" begin
        for row in rows
            @testset "$(row.check)" begin
                @test row.pass
            end
        end
    end
    println("Verification files: ", OUTPUT)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

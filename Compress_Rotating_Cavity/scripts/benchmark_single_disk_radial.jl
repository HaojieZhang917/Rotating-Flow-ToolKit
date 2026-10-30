#!/usr/bin/env julia

using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "SingleDiskRadial.jl"))
using .SingleDiskRadial

function csv_cell(value)
    value_string = string(value)
    return occursin(r"[,\"\n\r]", value_string) ?
           "\"" * replace(value_string, "\"" => "\"\"") * "\"" : value_string
end

function options(args)
    values = Dict("degree" => "80", "eta-max" => "20.0", "ds" => "0.02", "tolerance" => "2e-10",
                  "ratio" => "100.0", "output" => joinpath(ROOT, "benchmark_results", "single_disk_radial"))
    if "--help" in args
        println("julia --project=. scripts/benchmark_single_disk_radial.jl [--degree N] [--eta-max L] [--ds h] [--tolerance tol] [--ratio r/r0] [--output DIR]")
        return nothing
    end
    iseven(length(args)) || error("options require --name value pairs; see --help")
    for i in 1:2:length(args)
        startswith(args[i], "--") || error("expected an option, got $(args[i])")
        key = args[i][3:end]
        haskey(values, key) || error("unknown option: $(args[i])")
        values[key] = args[i+1]
    end
    ratio = parse(Float64, values["ratio"])
    ratio > 1 || error("radius ratio must exceed one")
    config = RadialConfig(degree=parse(Int, values["degree"]),
                          eta_max=parse(Float64, values["eta-max"]),
                          ds=parse(Float64, values["ds"]), s_end=log(ratio),
                          tolerance=parse(Float64, values["tolerance"]))
    return config, abspath(values["output"])
end

function decay_fit(run, floor)
    # Fixed, reported fit interval avoids both entry transients and numerical floor.
    selected = filter(m -> 0.8 <= m.s <= 1.8 && m.error_fg > max(100floor, 1e-10), run.metrics)
    length(selected) >= 4 || return NaN
    s = getproperty.(selected, :s)
    slope = hcat(ones(length(s)), s) \ log.(getproperty.(selected, :error_fg))
    return slope[2]
end

function write_results(runs, output)
    mkpath(output)
    history_fields = (:station, :s, :r_ratio, :error_fg, :error_f, :error_g,
                      :H_top, :wall_Fp, :wall_Gp, :residual, :iterations, :min_F,
                      :continuity_balance)
    open(joinpath(output, "history.csv"), "w") do io
        println(io, "case," * join(string.(history_fields), ","))
        for run in runs, m in run.metrics
            println(io, string(run.inlet.case), ",", join([string(getproperty(m, f)) for f in history_fields], ","))
        end
    end
    open(joinpath(output, "profiles.csv"), "w") do io
        println(io, "case,station,s,r_ratio,eta,F,G,H,F_ref,G_ref,H_ref")
        for run in runs
            ratios = [1.0, 1.5, 2.0, 5.0, 10.0, exp(run.s[end])]
            indices = sort(unique([argmin(abs.(run.s .- log(r))) for r in ratios if r <= exp(run.s[end]) + 1e-12]))
            ref = run.reference
            for i in indices, j in eachindex(ref.eta)
                println(io, join((run.inlet.case, i-1, run.s[i], exp(run.s[i]), ref.eta[j],
                                  run.F[j,i], run.G[j,i], run.H[j,i], ref.F[j], ref.G[j], ref.H[j]), ","))
            end
        end
    end
    floor = maximum(m.error_fg for m in runs[1].metrics)
    open(joinpath(output, "summary.csv"), "w") do io
        println(io, "case,status,stations,initial_error_fg,final_error_fg,final_over_initial,fitted_kappa_s0.8_to_1.8,max_residual,max_continuity_imbalance")
        for run in runs
            first_error, last_error = run.metrics[1].error_fg, run.metrics[end].error_fg
            reduction = first_error > 1e-15 ? last_error / first_error : NaN
            residuals = [m.residual for m in run.metrics if isfinite(m.residual)]
            balances = [abs(m.continuity_balance) for m in run.metrics if isfinite(m.continuity_balance)]
            println(io, join((run.inlet.case, run.status, length(run.s), first_error, last_error,
                              reduction, decay_fit(run, floor),
                              isempty(residuals) ? NaN : maximum(residuals),
                              isempty(balances) ? NaN : maximum(balances)), ","))
        end
    end
    open(joinpath(output, "inlets.csv"), "w") do io
        println(io, "case,amplitude_f,amplitude_g,center,power,description")
        for run in runs
            i = run.inlet
            println(io, join(csv_cell.((i.case, i.af, i.ag, i.center, i.power, i.description)), ","))
        end
    end
    ref = runs[1].reference
    open(joinpath(output, "reference.txt"), "w") do io
        println(io, "Single disk, quiescent exterior, no wall suction, zero external radial pressure gradient")
        println(io, "eta=z*sqrt(Omega/nu); s=log(r/r0)")
        println(io, "u_r=Omega*r*F; u_theta=Omega*r*G; u_z=sqrt(nu*Omega)*H")
        println(io, "degree=$(ref.config.degree), eta_max=$(ref.config.eta_max)")
        println(io, "method=$(runs[1].config.method), actual_ds=$(runs[1].ds)")
        println(io, "newton_tolerance=$(runs[1].config.tolerance), reverse_tolerance=$(runs[1].config.reverse_tolerance)")
        println(io, "reference_residual=$(ref.residual)")
        println(io, "F_prime_wall=$((ref.D*ref.F)[1])")
        println(io, "G_prime_wall=$((ref.D*ref.G)[1])")
        println(io, "H_top=$(ref.H[end])")
        println(io, "discrete_reference_max_error=$floor")
        for run in runs
            println(io, "$(run.inlet.case): $(run.status); $(run.message)")
        end
    end
    open(joinpath(output, "results.md"), "w") do io
        println(io, "# Single-disk radial attraction experiment\n")
        println(io, "Steady axisymmetric boundary-layer model; all cases retain the radial non-similar derivatives.\n")
        println(io, "Configuration: degree=$(ref.config.degree), eta_max=$(ref.config.eta_max), actual ds=$(runs[1].ds), requested r/r0=$(exp(runs[1].config.s_end)), Newton residual tolerance=$(runs[1].config.tolerance).\n")
        println(io, "| Inlet | Status | Initial relative FG error | Final relative FG error | Fitted spatial exponent |")
        println(io, "|---|---|---:|---:|---:|")
        for run in runs
            @printf(io, "| %s | %s | %.6e | %.6e | %.5f |\n", string(run.inlet.case),
                    string(run.status), run.metrics[1].error_fg, run.metrics[end].error_fg, decay_fit(run, floor))
        end
        println(io, "\nFitted exponents use 0.8 <= log(r/r0) <= 1.8; they are finite-interval fits, not eigenvalues. Values below max(100 times the control drift, 1e-10) are excluded.")
        println(io, "Late-stage error plateaus can reflect the nonlinear stopping tolerance. Errors are relative to the reference on the same finite eta grid; they are not absolute errors relative to the infinite-domain solution.")
        println(io, "The inlet has prescribed F,G only. Its H entries are NaN because inlet H is not independently prescribed or inferred by imposing F_s=0.")
        println(io, "For s>0, H is solved with H_eta+2F+F_s=0. H at the far boundary is not prescribed.")
        println(io, "\nThese finite inlet experiments test attraction within the steady axisymmetric outward-flow boundary-layer model. They do not prove global attraction for arbitrary inlets, finite disk NS flows, or three-dimensional physical stability.")
        println(io, "The r/r0 interval is a mathematical continuation interval; no physical local Reynolds-number interval has been assigned.")
        println(io, "No reverse velocity clipping, FLARE, or projection onto a similarity profile is used. A converged reverse flow terminates the run explicitly.")
        println(io, "\nVerification results are written separately by scripts/verify_single_disk_radial.jl.")
    end
end

function main(args=ARGS)
    parsed = options(args)
    parsed === nothing && return
    config, output = parsed
    BLAS.set_num_threads(1)
    reference = reference_profile(config)
    @printf("Reference: Fp(0)=%.12f, Gp(0)=%.12f, Htop=%.12f, residual=%.3e\n",
            (reference.D*reference.F)[1], (reference.D*reference.G)[1], reference.H[end], reference.residual)
    runs = []
    for name in (:similarity, :small_shape, :small_opposite, :large_radial, :large_swirl)
        run = march(reference, inlet_profile(reference; case=name))
        push!(runs, run)
        @printf("%-16s %-18s E: %.4e -> %.4e, r/r0=%.3f\n", string(name),
                string(run.status), run.metrics[1].error_fg, run.metrics[end].error_fg, exp(run.s[end]))
        flush(stdout)
    end
    write_results(runs, output)
    println("Saved: ", output)
    all(run.status == :completed for run in runs) || error("one or more cases stopped; see reference.txt")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

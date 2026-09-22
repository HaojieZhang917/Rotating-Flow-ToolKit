#!/usr/bin/env julia

using LinearAlgebra, Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const AUDIT = joinpath(ROOT, "work", "results", "presentation_evidence_audit_20260908")
const FULLROOT = joinpath(ROOT, "work", "results", "corrected_energy_epsilon_fold_chain")
const MODEROOT = joinpath(AUDIT, "tangent")
const OUTERFILE = joinpath(ROOT, "work", "results", "corrected_energy_outer_limit", "outer_profile.csv")
const OUTER_RO = -0.470367416464038

function csvdict(path)
    lines = filter(!isempty, strip.(readlines(path)))
    names = split(lines[1], ',')
    [Dict(names .=> parse.(Float64, split(line, ','))) for line in lines[2:end]]
end

function interp(x, y, xnew)
    out = similar(xnew)
    for i in eachindex(xnew)
        q = xnew[i]
        j = clamp(searchsortedlast(x, q), 1, length(x)-1)
        a = (q-x[j])/(x[j+1]-x[j])
        out[i] = (1-a)*y[j] + a*y[j+1]
    end
    out
end

function curve(rows, xname, yname, grid)
    clean = filter(r -> isfinite(r[xname]) && isfinite(r[yname]), rows)
    x = [r[xname] for r in clean]
    y = [r[yname] for r in clean]
    order = sortperm(x)
    interp(x[order], y[order], grid)
end

relerr(a, b) = (norm(a-b)/max(norm(b), 1e-14), maximum(abs, a-b))
token(epsilon) = replace(@sprintf("%.3g", epsilon), "."=>"p")

function table_row(path, epsilon)
    rows = csvdict(path)
    only(filter(r -> isapprox(r["epsilon"], epsilon; atol=1e-14, rtol=0), rows))
end

function main()
    mkpath(AUDIT)
    grid = collect(range(0, 10; length=1001))
    outer = csvdict(OUTERFILE)
    targets = Dict(
        "G" => curve(outer, "Y", "g", grid),
        "T" => curve(outer, "Y", "theta", grid),
        "H_over_epsilon" => curve(outer, "Y", "h", grid),
    )

    open(joinpath(AUDIT, "baseflow_outer_profile_errors.csv"), "w") do io
        println(io, "N,epsilon,field,relative_l2,max_abs,Y_range")
        for N in (240, 280, 320), epsilon in (0.02, 0.01)
            run = N == 240 ? "convergence_N240" : N == 280 ? "production_N280" : "production_N320"
            file = joinpath(FULLROOT, run, "profiles", "profile_epsilon_$(token(epsilon)).csv")
            rows = csvdict(file)
            converted = [Dict(
                "Y" => epsilon*r["eta"],
                "G" => r["G"],
                "T" => r["T"],
                "H_over_epsilon" => r["H"]/epsilon,
            ) for r in rows]
            for field in ("G", "T", "H_over_epsilon")
                a = curve(converted, "Y", field, grid)
                e, m = relerr(a, targets[field])
                @printf(io, "%d,%.12e,%s,%.12e,%.12e,0:10\n", N, epsilon, field, e, m)
            end
        end
    end

    fields = ("vG_scaled", "vT_scaled", "vH_scaled", "vF_scaled")
    open(joinpath(AUDIT, "zero_mode_outer_collapse.csv"), "w") do io
        println(io, "N,epsilon,reference_epsilon,field,relative_l2,max_abs,Y_range")
        for N in (240, 280, 320)
            refrows = csvdict(joinpath(MODEROOT, "N$(N)", "mode_epsilon_0p01.csv"))
            for epsilon in (0.02, 0.015)
                rows = csvdict(joinpath(MODEROOT, "N$(N)", "mode_epsilon_$(token(epsilon)).csv"))
                for field in fields
                    a = curve(rows, "Y", field, grid)
                    b = curve(refrows, "Y", field, grid)
                    dot(a, b) < 0 && (a .*= -1)
                    e, m = relerr(a, b)
                    @printf(io, "%d,%.12e,%.12e,%s,%.12e,%.12e,0:10\n",
                        N, epsilon, 0.01, field, e, m)
                end
            end
        end
    end

    run = "production_N320"
    table = joinpath(FULLROOT, run, "epsilon_fold_table.csv")
    open(joinpath(AUDIT, "inner_fixed_eta_scaling.csv"), "w") do io
        println(io, "epsilon,eta,G_over_epsilon,TminusTw_over_epsilon,F_over_epsilon2,H_over_epsilon2,core_compatibility,Ro_minus_outer,Tw_minus_outer")
        for epsilon in (0.1, 0.05, 0.02, 0.01, 0.005)
            rows = csvdict(joinpath(FULLROOT, run, "profiles", "profile_epsilon_$(token(epsilon)).csv"))
            row = table_row(table, epsilon)
            s = (2-row["Ro_f"]-row["Ro_f"]^2)/2
            omega = s+row["Ro_f"]
            chiw = 2-row["Tw_f"]
            defect = chiw*s^2-omega^2
            twouter = 1.66158741971899
            for eta in (1.0, 2.0, 5.0)
                g = curve(rows, "eta", "G", [eta])[1]
                t = curve(rows, "eta", "T", [eta])[1]
                f = curve(rows, "eta", "F", [eta])[1]
                h = curve(rows, "eta", "H", [eta])[1]
                @printf(io, "%.12e,%.8f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                    epsilon, eta, g/epsilon, (t-row["Tw_f"])/epsilon,
                    f/epsilon^2, h/epsilon^2, defect,
                    row["Ro_f"]-OUTER_RO, row["Tw_f"]-twouter)
            end
        end
    end
end

main()

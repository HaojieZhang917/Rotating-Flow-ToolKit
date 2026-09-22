#!/usr/bin/env julia

"""
Audit the radial transport direction in the existing generic-tail
first-order inner profiles.  This is a marching-admissibility diagnostic,
not a new base-flow solve.
"""

using DelimitedFiles
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const INPUT = joinpath(ROOT, "work", "results",
                       "generic_inner_first_order_20260909")
const OUTPUT = joinpath(ROOT, "work", "results",
                        "nonsimilar_heated_bek_20260909")

function read_profile(path)
    raw, header = readdlm(path, ',', Float64; header=true)
    names = vec(String.(header))
    columns = Dict(name => raw[:, i] for (i, name) in enumerate(names))
    columns
end

function zero_crossings(eta, velocity; eta_max=8.0, reltol=1e-6)
    keep = findall((eta .> 0) .& (eta .<= eta_max))
    e = eta[keep]
    u = velocity[keep]
    scale = maximum(abs, u)
    threshold = reltol * max(scale, eps(Float64))
    crossings = Float64[]
    last_index = nothing
    for j in eachindex(u)
        abs(u[j]) <= threshold && continue
        if last_index !== nothing && signbit(u[j]) != signbit(u[last_index])
            k = last_index
            push!(crossings, e[k] - u[k] * (e[j] - e[k]) / (u[j] - u[k]))
        end
        last_index = j
    end
    crossings, minimum(u), maximum(u), threshold
end

function main()
    mkpath(OUTPUT)
    cases = [
        (-0.50, "profile_Ro_0p500000.csv"),
        (-0.470367, "profile_Ro_0p470367.csv"),
        (-0.40, "profile_Ro_0p400000.csv"),
    ]
    rows = String[
        "Ro,source,eta_max,min_physical_coefficient,max_physical_coefficient," *
        "relative_zero_threshold,leading_amplitude,status," *
        "interior_zero_crossing_count,zero_crossings_eta"
    ]
    for (Ro, filename) in cases
        path = joinpath(INPUT, filename)
        data = read_profile(path)
        eta = data["eta"]
        # u_R/(Omega*R) = -Ro*F = epsilon*(-Ro*F1)+...
        velocity = -Ro .* data["F1"]
        crossings, umin, umax, threshold =
            zero_crossings(eta, velocity; eta_max=8.0)
        amplitude = max(abs(umin), abs(umax))
        status = amplitude < 1e-8 ? "first_order_degenerate" :
                 isempty(crossings) ? "one_sign_interior" :
                                      "sign_change_in_eta"
        status == "first_order_degenerate" && empty!(crossings)
        crossing_text = isempty(crossings) ? "" :
            join((@sprintf("%.8g", x) for x in crossings), ";")
        push!(rows, @sprintf("%.9g,%s,8,%.12e,%.12e,%.12e,%.12e,%s,%d,%s",
                            Ro, filename, umin, umax, threshold, amplitude,
                            status, length(crossings), crossing_text))
        @printf("Ro=%+.7f: min(-Ro F1)=%+.4e max=%+.4e status=%s crossings=%s\n",
                Ro, umin, umax, status,
                isempty(crossings) ? "none" : crossing_text)
    end
    open(joinpath(OUTPUT, "radial_direction_audit.csv"), "w") do io
        for row in rows
            println(io, row)
        end
    end
end

main()

#!/usr/bin/env julia

using Printf
using Plots

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const INPUT_ROOT = joinpath(ROOT, "work", "results",
    "presentation_evidence_audit_20260908", "rossby_dense")
const OUTPUT = joinpath(INPUT_ROOT, "summary")

function read_rows(path)
    lines = readlines(path)
    header = split(first(lines), ',')
    [Dict(header .=> split(line, ',')) for line in lines[2:end]]
end

function collect_stage(stage, degree)
    path = joinpath(INPUT_ROOT, "$(stage)_N$(degree)", "topology_table.csv")
    isfile(path) || return Dict{Tuple{String,Float64},String}()
    Dict((row["model"], parse(Float64, row["Ro"])) => row["topology"]
         for row in read_rows(path))
end

function main()
    mkpath(OUTPUT)
    sources = [
        (120, "stage1"), (120, "stage2"),
        (160, "stage3"), (200, "stage4c"),
        (240, "stage5"), (280, "stage6")]
    tables = Dict((n, stage) => collect_stage(stage, n) for (n, stage) in sources)
    rossby = [-0.30, -0.2875, -0.25, -0.20, -0.175, -0.1625, -0.15]
    degrees = [120, 160, 200, 240, 280]

    open(joinpath(OUTPUT, "second_window_grid_sensitivity.csv"), "w") do io
        println(io, "Ro," * join(["N$(n)" for n in degrees], ','))
        for ro in rossby
            values = String[]
            for n in degrees
                candidates = [table for ((degree, _), table) in tables if degree == n]
                value = "not_sampled"
                for table in candidates
                    value = get(table, ("consistent", ro), value)
                end
                push!(values, value)
            end
            @printf(io, "%.4f,%s\n", ro, join(values, ','))
        end
    end

    code = Dict("tail" => 0.0, "fold" => 1.0, "smooth" => 2.0)
    markers = Dict(120 => :circle, 160 => :rect, 200 => :diamond,
                   240 => :utriangle, 280 => :star5)
    plot(; xlabel="Rossby number", ylabel="classification code",
         yticks=([0,1,2], ["tail", "fold", "smooth/to Tw cap"]),
         xlims=(-0.31,-0.14), ylims=(-0.25,1.25), framestyle=:box,
         legend=:topright, gridalpha=0.25, size=(850,480))
    for n in degrees
        xs = Float64[]
        ys = Float64[]
        for ro in rossby
            value = "not_sampled"
            for ((degree, _), table) in tables
                degree == n || continue
                value = get(table, ("consistent", ro), value)
            end
            haskey(code, value) || continue
            push!(xs, ro); push!(ys, code[value])
        end
        scatter!(xs, ys; marker=markers[n], markersize=7, label="N=$n")
    end
    savefig(joinpath(OUTPUT, "second_window_grid_sensitivity.png"))
    println("Outputs: $OUTPUT")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

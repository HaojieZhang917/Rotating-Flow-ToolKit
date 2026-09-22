#!/usr/bin/env julia

"""Export representative corrected-energy branch data as Tecplot ASCII zones.

The source CSV files are continuation outputs from the dense Rossby audit:
* Ro = -1.0: traditional fold and consistent smooth branch;
* Ro = -0.5: traditional thermal tail and consistent fold.

Each continuation branch and pseudo-arclength segment is exported as a separate
zone. `segment_code=0` denotes fixed-Tw continuation and `segment_code=1`
denotes the pseudo-arclength segment that crosses a fold. Thermal length is
reported as NaN when absent from an arc source file.
"""

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SOURCE = joinpath(ROOT, "work", "results",
    "presentation_evidence_audit_20260908", "rossby_dense", "stage1_N120")
const OUTPUT = joinpath(ROOT, "work", "results",
    "presentation_evidence_audit_20260908", "tecplot")

const ZONES = [
    ("traditional_Ro_m1_fixed_Tw", "branch_traditional_Ro_m1p00.csv", 0),
    ("traditional_Ro_m1_pseudoarclength", "arc_traditional_Ro_m1p00.csv", 1),
    ("consistent_Ro_m1_fixed_Tw", "branch_consistent_Ro_m1p00.csv", 0),
    ("traditional_Ro_m0p5_fixed_Tw", "branch_traditional_Ro_m0p50.csv", 0),
    ("consistent_Ro_m0p5_fixed_Tw", "branch_consistent_Ro_m0p50.csv", 0),
    ("consistent_Ro_m0p5_pseudoarclength", "arc_consistent_Ro_m0p50.csv", 1),
]

function read_rows(path)
    lines = readlines(path)
    isempty(lines) && error("Empty source file: $path")
    header = split(first(lines), ',')
    rows = [Dict(header .=> split(line, ',')) for line in lines[2:end] if !isempty(line)]
    isempty(rows) && error("No data rows in source file: $path")
    return rows
end

value(row, name, fallback="NaN") = get(row, name, fallback)

function write_zone(io, title, rows, segment_code)
    println(io, "ZONE T=\"$(title)\", I=$(length(rows)), J=1, K=1, DATAPACKING=POINT")
    for row in rows
        println(io, join((
            value(row, "point"),
            value(row, "Ro"),
            value(row, "Tw"),
            value(row, "Hinf"),
            value(row, "residual"),
            value(row, "thermal_length"),
            string(segment_code),
        ), " "))
    end
end

function main()
    mkpath(OUTPUT)
    output = joinpath(OUTPUT, "representative_traditional_consistent_branches.dat")
    open(output, "w") do io
        println(io, "TITLE=\"Corrected-energy BEK representative continuation branches\"")
        println(io, "VARIABLES=\"point\" \"Ro\" \"Tw\" \"Hinf\" \"residual\" \"thermal_length\" \"segment_code\"")
        println(io, "# segment_code: 0=fixed-Tw continuation, 1=pseudo-arclength fold-crossing segment")
        println(io, "# NaN thermal_length: source arc files do not contain that diagnostic")
        for (title, filename, segment_code) in ZONES
            path = joinpath(SOURCE, filename)
            isfile(path) || error("Missing source data: $path")
            write_zone(io, title, read_rows(path), segment_code)
        end
    end
    println("Wrote: $output")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

#!/usr/bin/env julia

"""Plot the single-disk radial-marching benchmark using Julia and PyCall.

Run with `julia --project=. scripts/plot_single_disk_radial.jl [output_directory]`.
The numerical solver itself does not require Python or matplotlib. This optional
plotting script uses the project's existing PyCall dependency and a Python
installation providing matplotlib. No Python source file is generated.
"""

using PyCall

const ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_OUTPUT = joinpath(ROOT, "benchmark_results", "single_disk_radial")

struct SimpleCSV
    index::Dict{String, Int}
    rows::Vector{Vector{String}}
end

function read_simple_csv(path::AbstractString)
    isfile(path) || error("Missing benchmark data: $path. Run the radial benchmark first.")
    lines = filter(line -> !isempty(strip(line)), readlines(path))
    isempty(lines) && error("Empty CSV: $path")
    header = strip.(split(first(lines), ','))
    index = Dict(String(name) => i for (i, name) in enumerate(header))
    rows = Vector{String}[]
    for (line_number, line) in enumerate(lines[2:end])
        fields = String.(strip.(split(line, ',')))
        length(fields) == length(header) || error("Invalid CSV row $(line_number + 1) in $path")
        push!(rows, fields)
    end
    return SimpleCSV(index, rows)
end

function number_column(table::SimpleCSV, name::String, rows=eachindex(table.rows))
    haskey(table.index, name) || error("CSV is missing column $name")
    column = table.index[name]
    return [parse(Float64, table.rows[i][column]) for i in rows]
end

function case_rows(table::SimpleCSV, name::String)
    column = table.index["case"]
    return findall(row -> row[column] == name, table.rows)
end

function ordered_cases(table::SimpleCSV)
    present = unique(row[table.index["case"]] for row in table.rows)
    preferred = ["similarity", "small_shape", "small_opposite", "large_radial", "large_swirl"]
    return vcat(filter(name -> name in present, preferred),
                filter(name -> !(name in preferred), present))
end

const CASE_LABELS = Dict(
    "similarity" => "Similarity control",
    "small_shape" => "Small shape change",
    "small_opposite" => "Small opposite change",
    "large_radial" => "Large radial change",
    "large_swirl" => "Large swirl change",
)

const CASE_COLORS = Dict(
    "similarity" => "0.55",
    "small_shape" => "#0072B2",
    "small_opposite" => "#009E73",
    "large_radial" => "#D55E00",
    "large_swirl" => "#CC79A7",
)

function save_figure(fig, directory, basename)
    fig.savefig(joinpath(directory, basename * ".png"), dpi=220, bbox_inches="tight")
    fig.savefig(joinpath(directory, basename * ".pdf"), bbox_inches="tight")
end

function plot_error_decay(plt, history::SimpleCSV, directory)
    fig, ax = plt.subplots(figsize=(7.8, 5.2))
    for (i, name) in enumerate(ordered_cases(history))
        rows = case_rows(history, name)
        radius = number_column(history, "r_ratio", rows)
        error = number_column(history, "error_fg", rows)
        valid = findall(j -> isfinite(radius[j]) && radius[j] > 0 &&
                            isfinite(error[j]) && error[j] > 0, eachindex(radius))
        color = get(CASE_COLORS, name, "C$(mod(i - 1, 10))")
        label = get(CASE_LABELS, name, name)
        # Exact zeros have no log ordinate and are omitted rather than replaced
        # with invented data. The legend entry is retained for an all-zero control.
        ax.loglog(radius[valid], error[valid], color=color,
                  linestyle=name == "similarity" ? "--" : "-",
                  linewidth=name == "similarity" ? 1.3 : 1.8, label=label)
    end
    ax.set_xlabel("Radius ratio r / r0")
    ax.set_ylabel("Relative weighted L2 error in (F, G)")
    ax.set_title("Radial development of non-similar inlet profiles")
    ax.grid(true, which="both", linewidth=0.55, alpha=0.28)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2,
              frameon=false, fontsize=9)
    fig.subplots_adjust(bottom=0.25)
    fig.text(0.5, 0.015,
             "Steady axisymmetric boundary-layer model. Exact zero errors are omitted on the log scale.",
             ha="center", fontsize=8, color="0.35")
    save_figure(fig, directory, "error_decay")
    plt.close(fig)
end

function plot_profiles(plt, profiles::SimpleCSV, directory)
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 5.4), sharey=true)
    axf, axg, axh = axes
    cases = ordered_cases(profiles)
    isempty(cases) && error("Profiles CSV contains no cases")

    # The reference is stored at every snapshot; only one copy is needed.
    first_rows = case_rows(profiles, first(cases))
    station = number_column(profiles, "station", first_rows)
    reference_rows = first_rows[findall(==(minimum(station)), station)]
    eta = number_column(profiles, "eta", reference_rows)
    order = sortperm(eta)
    eta = eta[order]
    reference_rows = reference_rows[order]
    axf.plot(number_column(profiles, "F_ref", reference_rows), eta,
             color="black", linestyle=":", linewidth=2.1, label="von Karman reference")
    axg.plot(number_column(profiles, "G_ref", reference_rows), eta,
             color="black", linestyle=":", linewidth=2.1)
    axh.plot(number_column(profiles, "H_ref", reference_rows), eta,
             color="black", linestyle=":", linewidth=2.1)

    for (i, name) in enumerate(cases)
        rows = case_rows(profiles, name)
        stations = number_column(profiles, "station", rows)
        inlet_rows = rows[findall(==(minimum(stations)), stations)]
        final_rows = rows[findall(==(maximum(stations)), stations)]
        color = get(CASE_COLORS, name, "C$(mod(i - 1, 10))")
        label = get(CASE_LABELS, name, name)
        for (ax, field) in [(axf, "F"), (axg, "G")]
            ax.plot(number_column(profiles, field, inlet_rows),
                    number_column(profiles, "eta", inlet_rows), color=color,
                    linestyle="--", linewidth=1.2, alpha=0.8)
            ax.plot(number_column(profiles, field, final_rows),
                    number_column(profiles, "eta", final_rows), color=color,
                    linestyle="-", linewidth=1.6, label=field == "F" ? label : "_nolegend_")
        end
        # H is determined by marching plus continuity, not independently supplied
        # at the inlet. Never plot the inlet H guess or its intentionally NaN data.
        h = number_column(profiles, "H", final_rows)
        he = number_column(profiles, "eta", final_rows)
        valid = findall(j -> isfinite(h[j]) && isfinite(he[j]), eachindex(h))
        axh.plot(h[valid], he[valid], color=color, linewidth=1.6)
    end
    for (ax, field) in [(axf, "F"), (axg, "G"), (axh, "H")]
        ax.set_xlabel(field)
        ax.set_title(field == "H" ? "Axial velocity (final only)" :
                     field == "F" ? "Radial velocity" : "Azimuthal velocity")
        ax.grid(true, linewidth=0.55, alpha=0.28)
    end
    axf.set_ylabel("Wall-normal coordinate eta")
    # The full domain remains in the CSV. A near-wall view makes the different
    # inlet shapes legible instead of compressing the boundary layer to one edge.
    axf.set_ylim(0, min(8.0, maximum(eta)))
    handles, labels = axf.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.045),
               ncol=3, frameon=false, fontsize=9)
    fig.suptitle("Scaled profiles: dashed inlet, solid final computed station", fontsize=12)
    fig.subplots_adjust(bottom=0.23, top=0.87, wspace=0.22)
    fig.text(0.5, 0.01,
             "Black dotted: reference. Near-wall view; full eta domain is retained in profiles.csv.",
             ha="center", fontsize=8, color="0.35")
    save_figure(fig, directory, "profiles")
    plt.close(fig)
end

function main(args=ARGS)
    length(args) <= 1 || error("Usage: plot_single_disk_radial.jl [output_directory]")
    directory = isempty(args) ? DEFAULT_OUTPUT : abspath(first(args))
    history = read_simple_csv(joinpath(directory, "history.csv"))
    profiles = read_simple_csv(joinpath(directory, "profiles.csv"))
    matplotlib = pyimport("matplotlib")
    matplotlib.use("Agg")
    plt = pyimport("matplotlib.pyplot")
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["pdf.fonttype"] = 42
    plot_error_decay(plt, history, directory)
    plot_profiles(plt, profiles, directory)
    println("Saved error_decay.png/.pdf and profiles.png/.pdf to ", directory)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

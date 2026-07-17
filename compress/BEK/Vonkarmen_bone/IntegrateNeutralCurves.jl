using Printf

const ROOT = @__DIR__
const BATCH_DIR = joinpath(ROOT,"neutral_curve_batch")
const OUTPUT_DIR = joinpath(ROOT,"neutral_curve_integrated")

const LEGACY_SOURCES = (
    (
        "ome=0.0_Tw=1.0.dat",
        "ome=0.0_Tw=1.0_model=lopez.dat",
        "legacy Lopez curve; coarse neutral tolerance",
    ),
    (
        "ome=0.0_Tw=1.0_compressible_Mr=0.3.dat",
        "ome=0.0_Tw=1.0_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat",
        "legacy full-compressible curve; coarse neutral tolerance",
    ),
    (
        "ome=0.0_Tw=1.02_model=lopez_propPert=off_baseProp=variable.dat",
        "ome=0.0_Tw=1.02_model=lopez.dat",
        "renamed to the standard Lopez convention",
    ),
    (
        "ome=0.0_Tw=1.02_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat",
        "ome=0.0_Tw=1.02_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat",
        "legacy parent-directory source",
    ),
    (
        "ome=0.0_Tw=1.02_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=variable.dat",
        "ome=0.0_Tw=1.02_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=variable.dat",
        "legacy parent-directory source",
    ),
    (
        "ome=0.0_Tw=1.02_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=frozen.dat",
        "ome=0.0_Tw=1.02_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=frozen.dat",
        "legacy parent-directory source",
    ),
    (
        "ome=0.0_Tw=1.04_model=lopez_propPert=off_baseProp=variable.dat",
        "ome=0.0_Tw=1.04_model=lopez.dat",
        "renamed to the standard Lopez convention",
    ),
    (
        "ome=0.0_Tw=1.04_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat",
        "ome=0.0_Tw=1.04_model=compressible_Mr=0.3_" *
        "propPert=on_baseProp=variable.dat",
        "legacy parent-directory source",
    ),
    (
        "ome=0.0_Tw=1.04_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=variable.dat",
        "ome=0.0_Tw=1.04_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=variable.dat",
        "legacy parent-directory source",
    ),
    (
        "ome=0.0_Tw=1.04_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=frozen.dat",
        "ome=0.0_Tw=1.04_model=compressible_Mr=0.3_" *
        "propPert=off_baseProp=frozen.dat",
        "legacy parent-directory source",
    ),
    (
        "ome=0.0_Tw=1.05.dat",
        "ome=0.0_Tw=1.05_model=lopez.dat",
        "legacy Lopez-only curve; no paired compressible file",
    ),
)

function load_curve(path)
    rows = Vector{Vector{Float64}}()
    for line in Iterators.drop(eachline(path),2)
        isempty(strip(line)) && continue
        push!(rows,parse.(Float64,split(strip(line))))
    end
    isempty(rows) && return zeros(0,7)
    return reduce(vcat,permutedims.(rows))
end

function curve_quality(path)
    data = load_curve(path)
    points = size(data,1)
    points == 0 && return (
        points=0,R_min=NaN,R_max=NaN,beta_min=NaN,beta_max=NaN,
        max_residual=NaN,direction="none",quality="invalid",
    )
    beta_difference = diff(data[:,3])
    direction = all(beta_difference .> 0) ? "increasing" :
        all(beta_difference .< 0) ? "decreasing" : "mixed"
    max_residual = maximum(min.(abs.(data[:,5]),abs.(data[:,7])))
    is_branch = occursin("_branch=",basename(path))
    quality = if is_branch && points >= 10 && max_residual <= 1.0e-6 &&
                 direction != "mixed"
        "validated_branch"
    elseif points >= 50 && max_residual <= 1.0e-6 &&
           direction != "mixed" && maximum(data[:,3]) >= 0.08
        "validated_strict"
    elseif points >= 50 && max_residual <= 5.0e-3 && direction != "mixed"
        "legacy_coarse"
    else
        "invalid"
    end
    return (
        points=points,R_min=minimum(data[:,2]),R_max=maximum(data[:,2]),
        beta_min=minimum(data[:,3]),beta_max=maximum(data[:,3]),
        max_residual=max_residual,direction=direction,quality=quality,
    )
end

function extract_temperature(name)
    matched = match(r"Tw=([0-9.]+)",name)
    matched === nothing && return NaN
    return parse(Float64,matched.captures[1])
end

mkpath(OUTPUT_DIR)
sources = Dict{String,Tuple{String,String}}()

for source in sort(filter(
    path -> begin
        name = basename(path)
        isfile(path) && endswith(name,".dat") &&
            !endswith(name,"_allbranches.dat") &&
            !startswith(name,"lopez_allbranches_")
    end,
    readdir(BATCH_DIR; join=true),
))
    destination = joinpath(OUTPUT_DIR,basename(source))
    cp(source,destination; force=true)
    sources[basename(destination)] = (source,"validated neutral_curve_batch source")
end

for (source_name,destination_name,note) in LEGACY_SOURCES
    source = joinpath(ROOT,source_name)
    destination = joinpath(OUTPUT_DIR,destination_name)
    if isfile(source)
        cp(source,destination; force=true)
        sources[destination_name] = (source,note)
    elseif isfile(destination)
        sources[destination_name] = (
            destination,note * "; preserved existing integrated copy",
        )
    else
        error(
            "Missing both legacy source and integrated copy: " *
            "$source, $destination",
        )
    end
end

ordered_files = sort(
    collect(keys(sources)); by=name -> (extract_temperature(name),name),
)

combined_path = joinpath(OUTPUT_DIR,"neutral_curves_all.dat")
open(combined_path,"w") do io
    println(
        io,
        "Variables=\"omega\" \"R\" \"beta\" \"alpha_r_1\" " *
        "\"alpha_i_1\" \"alpha_r_2\" \"alpha_i_2\"",
    )
    for name in ordered_files
        println(io,"Zone T=\"$(chop(name; tail=4))\"")
        for line in Iterators.drop(eachline(joinpath(OUTPUT_DIR,name)),2)
            isempty(strip(line)) || println(io,line)
        end
    end
end

manifest_path = joinpath(OUTPUT_DIR,"manifest.tsv")
open(manifest_path,"w") do io
    println(
        io,
        "file\tTw\tquality\tpoints\tR_min\tR_max\tbeta_min\tbeta_max\t" *
        "max_neutral_residual\tbeta_direction\tsource\tnote",
    )
    for name in ordered_files
        source,note = sources[name]
        quality = curve_quality(joinpath(OUTPUT_DIR,name))
        fields = (
            name,extract_temperature(name),quality.quality,quality.points,
            quality.R_min,quality.R_max,quality.beta_min,quality.beta_max,
            quality.max_residual,quality.direction,relpath(source,ROOT),note,
        )
        println(io,join(fields,'\t'))
    end
end

println("Integrated $(length(sources)) curve files into $OUTPUT_DIR")
println("Manifest: $manifest_path")
println("Tecplot multi-zone file: $combined_path")

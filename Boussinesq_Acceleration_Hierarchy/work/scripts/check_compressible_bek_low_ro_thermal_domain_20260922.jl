using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "work", "src", "CompressibleBEKLinearStability.jl"))
const CBEK = CompressibleBEKLinearStability
const OUTDIR = joinpath(ROOT,"work","results",
    "compressible_bek_thermal_limit_gate_20260922")
mkpath(OUTDIR)

function value_at(source, field, y)
    index = round(Int, y / (source.y[2]-source.y[1])) + 1
    return getfield(source,field)[index]
end

rows = NamedTuple[]
for (Ro,lengths) in ((-0.2,(40.0,80.0)),
                     (-0.1,(40.0,80.0,160.0)),
                     (-0.05,(40.0,80.0,160.0,320.0)))
    Mr = abs(Ro)*0.3
    for L in lengths
        source = CBEK.source_profiles(Ro,0.72;ymax=L)
        fscale = maximum(abs,source.f)
        source_tail = abs(value_at(source,:f,0.75L))/fscale
        row = (Ro=Ro,Mr=Mr,L=L,
            decay_exponent=0.72*Ro*last(source.w),
            decay_products=-0.72*Ro*last(source.w)*L,
            f_y1=value_at(source,:f,1.0),
            f_y10=value_at(source,:f,10.0),
            f_y20=value_at(source,:f,20.0),
            f_tail_fraction=source_tail,
            max_temperature_change=0.2*Mr^2*fscale,
            q_y10=value_at(source,:q,10.0),
            q_tail=value_at(source,:q,0.75L))
        push!(rows,row)
        @printf("Ro=%6.3f L=%5.0f decay_lengths=%.3f f(1)=%.6f f_tail/max=%.3e Tdev=%.3e\n",
            Ro,L,row.decay_products,row.f_y1,row.f_tail_fraction,
            row.max_temperature_change)
        flush(stdout)
    end
end
open(joinpath(OUTDIR,"low_ro_thermal_domain.csv"),"w") do io
    println(io,join(keys(first(rows)),','))
    foreach(r->println(io,join(values(r),',')),rows)
end

# Two example profiles for the gap-filter figure: evanescent (Type-I-like) and inertial (critical level) cases.
src = read(joinpath(@__DIR__, "analyze_gap_filter.jl"), String)
eval(Meta.parseall(src[1:findfirst("Re0, R0 = 1000.0", src)[1]-1]))
b0 = cavity_baseflow(1000.0; N=300); R0 = 330.0; K0 = bary_interp(b0.z, b0.G, 0.5); zr0, zs0 = edges(b0)
open(joinpath(@__DIR__, "results", "gap_filter", "profiles.csv"), "w") do fh
    println(fh, "case,z,absu,G,sigma_local_over_f")
    for (name, a, n, om) in (("evanescent_a0.45_n16_wbar8", 0.45, 16, 8/R0), ("inertial_a0.2_n-10_c0.7", 0.2, -10, 0.7*10/R0))
        beta = n/R0; c = -om/beta
        op, q = forced(b0, R0, a, n, c)
        for t in range(0.0, 1.0, length=401)
            G = bary_interp(b0.z, b0.G, t)
            println(fh, join((name, t, abs(bary_interp(op.z, comp(op, q, 1), t)), G, (om + beta*(1 - G))/(2G/R0)), ","))
        end
    end
end
println("zr=$zr0 zs=$zs0 K=$K0")

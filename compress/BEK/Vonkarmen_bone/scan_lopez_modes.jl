include(joinpath(@__DIR__, "LopezStability.jl"))
include(joinpath(@__DIR__, "LopezBaseflow.jl"))
include(joinpath(@__DIR__, "CRD_STA.jl"))

using .LopezStability
using .LopezBaseflow
using .CRD_BF
using BSplineKit

Tw = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.05
R = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 402.0
beta = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0481
N_cheb = 99

profile = solve_lopez_baseflow(Tw)
D, D2, z_cheb = CRD_BF.Cheb(N_cheb)
sample(values) = begin
    spline = BSplineKit.interpolate(profile.z, values, BSplineOrder(4))
    vec([spline(point) for point in z_cheb])
end
F, G, H, T = sample(profile.F), sample(profile.G), sample(profile.H), sample(profile.T)

found = NamedTuple[]
for shift in 0.025:0.025:0.800
    try
        values, vectors, info = eigsol_lopez(
            F, G, H, T, R, 0.0, beta, N_cheb, D, D2,
            ComplexF64(shift), 2;
            tol=1.0e-10,
            maxit=1000,
            return_info=true,
        )
        for index in eachindex(values)
            value = values[index]
            residual = info.residuals[index]
            if isfinite(real(value)) && isfinite(imag(value)) &&
               -0.1 <= real(value) <= 1.0 && abs(imag(value)) <= 0.15 &&
               residual <= 1.0e-8 &&
               all(abs(value - mode.alpha) > 1.0e-5 for mode in found)
                push!(found, (
                    alpha=value, residual=residual,
                    vector=copy(vectors[:, index]), shift=shift,
                ))
            end
        end
    catch exception
        println("shift=$shift failed: $(sprint(showerror, exception))")
    end
end

sort!(found; by=mode -> abs(imag(mode.alpha)))
println("Tw=$Tw R=$R beta=$beta: $(length(found)) unique near-real modes")
for (index, mode) in enumerate(found)
    println(
        "mode=$index alpha=$(mode.alpha) residual=$(mode.residual) " *
        "source_shift=$(mode.shift)",
    )
end

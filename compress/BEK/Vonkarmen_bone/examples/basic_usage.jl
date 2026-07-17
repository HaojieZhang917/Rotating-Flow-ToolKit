using Pkg

Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using RotatingDiskFlow
using Plots

profile = RotatingDiskFlow.LopezBaseflow.solve_lopez_baseflow(1.10)

plot(
    profile.F,
    profile.z;
    xlabel="F, G, H",
    ylabel="z",
    label="F",
    linewidth=2,
)
plot!(profile.G, profile.z; label="G", linewidth=2)
plot!(profile.H, profile.z; label="H", linewidth=2)

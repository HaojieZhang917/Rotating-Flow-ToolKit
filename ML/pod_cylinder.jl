# Run from this folder:
#   julia --project=. -e 'using Pkg; Pkg.instantiate()'
#   julia --project=. pod_cylinder.jl
# The data come from the original PINNs repository, not a neural-network prediction.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")  # save figures without opening windows
using MAT, LinearAlgebra, Statistics, Printf, Downloads, Plots

# Change these two settings to explore the decomposition yourself.
const RETAINED_MODES = 2
const REMOVE_MEAN = true

const BASE = @__DIR__
const OUT = joinpath(BASE, "results")
const DATA = joinpath(BASE, "cylinder_nektar_wake.mat")
# const DATA_URL = "https://raw.githubusercontent.com/maziarraissi/PINNs/master/main/Data/cylinder_nektar_wake.mat"
isfile(DATA) || Downloads.download(DATA_URL, DATA)
mkpath(OUT)

# U_star: spatial point x velocity component x time.
d = matread(DATA)
U = d["U_star"]
xy = d["X_star"]
t = vec(d["t"])
n, ncomponents, nt = size(U)
@assert ncomponents == 2 && size(xy) == (n, 2) && length(t) == nt
@assert all(isfinite, U)

# The first n rows are streamwise velocity; the next n are transverse velocity.
X = vcat(U[:, 1, :], U[:, 2, :])
mu = REMOVE_MEAN ? mean(X; dims=2) : zeros(size(X, 1), 1)
Xc = X .- mu
F = svd(Xc; full=false)
energy = F.S .^ 2
sum(energy) > 0 || error("The selected data have no fluctuation energy.")
ratio = energy ./ sum(energy)
cumulative = cumsum(ratio)
r = RETAINED_MODES
1 <= r <= length(F.S) || error("RETAINED_MODES is out of range.")
Phi = F.U[:, 1:r]
A = Diagonal(F.S[1:r]) * F.Vt[1:r, :]
Xrec = mu .+ Phi * A
relative_error = norm(Xc - Phi * A) / norm(Xc)

# Independent consistency checks: orthogonality and the truncated-SVD error identity.
@assert isapprox(Phi' * Phi, Matrix{Float64}(I, r, r); atol=1e-10)
@assert isapprox(relative_error^2, sum(ratio[r+1:end]); atol=1e-10)

println("Data: ", n, " spatial points, ", nt, " snapshots; X size = ", size(X))
println("The time coordinate is nondimensional, not seconds.")
for k in 1:min(10, length(ratio))
    @printf("Mode %2d: %9.5f%%; cumulative: %9.5f%%\n", k, 100ratio[k], 100cumulative[k])
end
@printf("Retained %d modes: %.5f%% captured, %.7f relative fluctuation error\n",
    r, 100cumulative[r], relative_error)

# Store the first 20 modes so they can be inspected and reconstructed later.
nsave = min(20, length(F.S))
matwrite(joinpath(OUT, "pod_results.mat"), Dict(
    "xy" => xy, "t" => t, "mean_field" => mu,
    "Phi" => F.U[:, 1:nsave],
    "A" => Diagonal(F.S[1:nsave]) * F.Vt[1:nsave, :],
    "energy_ratio" => ratio, "singular_values" => F.S,
    "retained_modes" => r, "relative_fluctuation_error" => relative_error
))

# Use the coordinates to place values on a grid; do not assume a flattening order.
xs, ys = sort(unique(xy[:, 1])), sort(unique(xy[:, 2]))
@assert length(xs) * length(ys) == n
@assert isapprox(diff(xs), fill(xs[2] - xs[1], length(xs) - 1); atol=1e-10)
@assert isapprox(diff(ys), fill(ys[2] - ys[1], length(ys) - 1); atol=1e-10)
ix, iy = Dict(v => j for (j, v) in enumerate(xs)), Dict(v => j for (j, v) in enumerate(ys))
function spatial_grid(values)
    z = fill(NaN, length(ys), length(xs))
    for j in 1:n
        z[iy[xy[j, 2]], ix[xy[j, 1]]] = values[j]
    end
    return z
end
function fieldplot(values, title_text; colorlimit=maximum(abs, values))
    return heatmap(xs, ys, spatial_grid(values); color=:RdBu,
        clims=(-colorlimit, colorlimit), aspect_ratio=:equal,
        xlims=(first(xs), last(xs)), ylims=(first(ys), last(ys)), widen=false,
        xlabel="x/D", ylabel="y/D", title=title_text, framestyle=:box)
end
default(fontfamily="sans-serif", dpi=150, titlefontsize=11, guidefontsize=10,
    tickfontsize=9, legendfontsize=9, left_margin=6Plots.mm,
    right_margin=4Plots.mm, bottom_margin=7Plots.mm, top_margin=4Plots.mm)

nk = min(12, length(ratio))
p1 = bar(1:nk, 100 .* ratio[1:nk]; xlabel="Mode", ylabel="Energy (%)",
    title="Individual energy", legend=false)
p2 = plot(1:nk, 100 .* cumulative[1:nk]; marker=:circle,
    xlabel="Modes retained", ylabel="Cumulative energy (%)", ylim=(0, 102),
    title="Cumulative energy", legend=false)
savefig(plot(p1, p2; size=(1000, 380)), joinpath(OUT, "01_energy.png"))

panels = [fieldplot(F.U[n+1:2n, k], @sprintf("Mode %d: %.2f%% (v component)", k, 100ratio[k])) for k in 1:4]
savefig(plot(panels...; layout=(2, 2), size=(1100, 650)), joinpath(OUT, "02_spatial_modes.png"))

A2 = Diagonal(F.S[1:2]) * F.Vt[1:2, :]
p1 = plot(t, A2[1, :]; label="a1", xlabel="t (nondimensional)", ylabel="Coefficient",
    title="Leading temporal coefficients")
plot!(p1, t, A2[2, :]; label="a2")
p2 = plot(A2[1, :], A2[2, :]; xlabel="a1", ylabel="a2", legend=false,
    title="Leading-pair phase portrait", aspect_ratio=:equal)
savefig(plot(p1, p2; size=(1100, 470)), joinpath(OUT, "03_time_and_phase.png"))

j = nt ÷ 2
original, reconstructed = X[n+1:2n, j], Xrec[n+1:2n, j]
limit = max(maximum(abs, original), maximum(abs, reconstructed))
p1 = fieldplot(original, "Original transverse velocity"; colorlimit=limit)
p2 = fieldplot(reconstructed, "Reconstruction with $r modes"; colorlimit=limit)
p3 = fieldplot(original - reconstructed, "Residual (separate color scale)")
savefig(plot(p1, p2, p3; layout=(3, 1), size=(1000, 1300)), joinpath(OUT, "04_reconstruction.png"))
println("Saved figures and pod_results.mat to: ", OUT)

#=
compare_turbulence_closures.jl
==============================

Closure-independence check (step (a)) of the turbulent-stator base flow, Re_h = 1000, local similarity:
  ML   : Michel mixing length + van Driest + Klebanoff (cavity_baseflow_turb default, edge = :peak)
  ML99 : same, outer scale = 0.99K thickness (edge = :first99)
  SA   : Spalart–Allmaras (noft2, boundary-layer shear), production off in the rotor half, nu_t confined
         to the stator layer by the same Klebanoff factor (rotor layer laminar)
  SAu  : SA without confinement (nu~ advected across the core into the rotor layer) - shows the
         opposite K shift when the rotor layer is (spuriously) turbulent.
R = 100 ... 600 (Omega r^2/nu = R^2 = 1e4 ... 3.6e5); SA continued in R (previous chi as seed).
The momentum-integral check (laminar Ekman + 1/7-law turbulent Bödewadt) was dropped: with the
Kármán profile family u = -eps v_c eta^{1/7}(1-eta), v = v_c eta^{1/7} the radial integral requires
eps^2 (0.504 - 0.746) = 2/9, i.e. it has no solution.
Output: results/turbulent_baseflow/closures.{csv,md}, closures.png (plot_turbulence_closures.py)
=#
include(joinpath(@__DIR__, "CavityLST.jl")); using .CavityLST
using LinearAlgebra, Printf
BLAS.set_num_threads(min(8, Sys.CPU_THREADS))
const OUT = joinpath(@__DIR__, "results", "turbulent_baseflow"); mkpath(OUT)
const Re_h, N = 1000.0, 200
Rs = collect(100.0:50.0:600.0)
rows = String[]; profs = Dict{String,Any}()
chi = nothing; chiu = nothing
for R in Rs
    global chi, chiu
    ml = cavity_baseflow_turb(Re_h, R; N=N)
    m99 = cavity_baseflow_turb(Re_h, R; N=N, edge=:first99)
    sa = cavity_baseflow_turb(Re_h, R; N=N, model=:sa, chi0=chi); chi = sa.chi
    sau = try
        r = cavity_baseflow_turb(Re_h, R; N=N, model=:sa, confine=false, chi0=chiu); chiu = r.chi; r
    catch
        nothing
    end
    for (name, r) in (("ML", ml), ("ML99", m99), ("SA", sa), ("SAu", sau))
        r === nothing && continue
        push!(rows, join((name, R, r.K, r.utau, r.delta_s_peak, maximum(r.nut), r.base.Pi), ","))
        R == 330.0 && (profs[name] = r)
    end
    @printf("R=%.0f  K: ML %.4f  ML99 %.4f  SA %.4f  SAu %s   u_tau: ML %.4f SA %.4f   max nu_t: ML %.2f SA %.2f\n", R, ml.K, m99.K, sa.K,
            sau === nothing ? "fail" : @sprintf("%.4f", sau.K), ml.utau, sa.utau, maximum(ml.nut), maximum(sa.nut))
end
open(joinpath(OUT, "closures.csv"), "w") do f
    println(f, "model,R,K,utau,delta_s_peak,max_nut,Pi"); foreach(r -> println(f, r), rows)
end
lam = cavity_baseflow(Re_h; N=N)
open(joinpath(OUT, "closures_profiles.csv"), "w") do f
    println(f, "model,z,G,F,nut")
    for j in eachindex(lam.z); println(f, "laminar,$(lam.z[j]),$(lam.G[j]),$(lam.F[j]),0"); end
    for (name, r) in profs, j in eachindex(r.base.z)
        println(f, "$name,$(r.base.z[j]),$(r.base.G[j]),$(r.base.F[j]),$(r.nut[j])")
    end
end
println("wrote ", OUT)

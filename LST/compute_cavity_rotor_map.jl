#=
compute_cavity_rotor_map.jl
===========================
Spatial growth-rate map of ROTOR-layer modes in the rotor–stator cavity (CavityLST, laminar Batchelor,
Re_h = 1000, rotor frame), over (n, omega_bar = omega R) at fixed R = 330 and 150.
All candidate eigenvalues in the window 0.05 < Re alpha < 1.2, -0.3 < Im alpha < 0.05 whose |v| peaks in
the rotor half (z < 0.25) are stored (first attempt kept only the most negative Im alpha and was swamped by
spurious |alpha| ~ 3 branch modes with Im alpha ~ -3). N = 89, P_N-P_{N-2} (Type-I target identical to N = 199).
Parallel: Distributed workers (BLAS threads contend within one process).
Output: results/cavity_rotor_map/candidates_R<R>.csv
=#
using Distributed
nw = min(16, Sys.CPU_THREADS); nprocs() == 1 && addprocs(nw)
@everywhere include(joinpath($(@__DIR__), "CavityLST.jl"))
@everywhere using .CavityLST
@everywhere using LinearAlgebra
@everywhere BLAS.set_num_threads(1)
@everywhere begin
    const B = cavity_baseflow(1000.0; N=300)
    function point(R, n, wb)
        op = cavity_pencil(B, R, n/R, wb/R; N=89, pressure=:PN2)
        out = Tuple{Float64,Float64,Float64}[]
        for m in cavity_spectrum(op; shift=0.4 + 0im, maxabs=3.0)
            a = m.alpha
            (0.05 < real(a) < 1.2 && -0.3 < imag(a) < 0.05) || continue
            x = full_state(op, m.vector); v = abs.(x[op.n+1:2op.n]); zp = op.z[argmax(v)]
            zp < 0.25 && push!(out, (real(a), imag(a), zp))
        end
        return (R, n, wb, out)
    end
end
const OUT = joinpath(@__DIR__, "results", "cavity_rotor_map"); mkpath(OUT)
for R in (330.0, 150.0)
    t0 = time()
    jobs = [(R, n, wb) for n in -30:2:60 for wb in -10:1:16]
    res = pmap(j -> point(j...), jobs)
    open(joinpath(OUT, "candidates_R$(Int(R)).csv"), "w") do f
        println(f, "R,n,wbar,alpha_r,alpha_i,z_peak")
        for (R_, n, wb, out) in res, (ar, ai, zp) in out
            println(f, join((R_, n, wb, ar, ai, zp), ","))
        end
    end
    println("R=$(R) done in $(round(time() - t0)) s"); flush(stdout)
end

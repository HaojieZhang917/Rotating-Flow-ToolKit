# Compressible rotating-disk base-flow marching solver with Sutherland viscosity.
#
# The solver uses the steady axisymmetric compressible boundary-layer equations
# in a rotating frame. Velocities are written as
#
#     u(r,z) = r * U(r,z),   v(r,z) = r * V(r,z),   w(r,z) = W(r,z),
#
# where V -> -1 in the far field. The model keeps ideal-gas density,
# Sutherland viscosity, variable conductivity, radial non-similarity terms, and
# viscous dissipation in the energy equation.
#
# Unknowns at each radius are U,V,W,T on a uniform physical z grid. Each radial
# station is solved implicitly with Newton iteration and finite-difference
# Jacobians, so the file only depends on Julia standard libraries.

module SutherlandMarching

using LinearAlgebra
using DelimitedFiles
using Printf
using Dates

export Params, Station, solve_marching, write_outputs, suth_mu, suth_k, suth_rho, main

Base.@kwdef struct Params
    Tw::Float64 = 1.5
    sigma::Float64 = 0.72
    gamma::Float64 = 1.4
    Minf::Float64 = 0.003       # local Mach is Mr(r)=Minf*r
    Tinf_dim::Float64 = 273.0
    S_dim::Float64 = 114.0
    r0::Float64 = 1.0
    rmax::Float64 = 100.0
    Nr::Int = 41
    zmax::Float64 = 20.0
    Nz::Int = 81
    newton_tol::Float64 = 1.0e-8
    newton_maxiter::Int = 14
    fd_eps::Float64 = 1.0e-6
    out_dir::String = "sutherland_marching_data"
    verbose::Bool = true
end

struct Station
    r::Float64
    U::Vector{Float64}
    V::Vector{Float64}
    W::Vector{Float64}
    T::Vector{Float64}
end

function suth_mu(T, p::Params)
    s = p.S_dim / p.Tinf_dim
    Ts = max(T, 1.0e-10)
    return Ts^(1.5) * (1.0 + s) / (Ts + s)
end

function suth_k(T, p::Params)
    return suth_mu(T, p) / p.sigma
end

suth_rho(T) = 1.0 / max(T, 1.0e-10)

function pack(U, V, W, T)
    return vcat(U, V, W, T)
end

function unpack(x::Vector{Float64}, Nz::Int)
    U = @view x[1:Nz]
    V = @view x[Nz+1:2Nz]
    W = @view x[2Nz+1:3Nz]
    T = @view x[3Nz+1:4Nz]
    return U, V, W, T
end

function initial_guess(z::Vector{Float64}, p::Params)
    U = 0.51023 .* z .* exp.(-0.8 .* z)
    V = -1.0 .+ exp.(-0.6159 .* z)
    W = -0.88447 .* (1.0 .- exp.(-0.8 .* z))
    T = 1.0 .+ (p.Tw - 1.0) .* exp.(-0.5 .* z)
    U[1] = 0.0
    V[1] = 0.0
    W[1] = 0.0
    T[1] = p.Tw
    U[end] = 0.0
    V[end] = -1.0
    T[end] = 1.0
    return pack(U, V, W, T)
end

function derivative_center(f, j, dz)
    return (f[j + 1] - f[j - 1]) / (2.0 * dz)
end

function diffusive_flux_second(q, coeff, j, dz)
    coeff_p = 0.5 * (coeff[j] + coeff[j + 1])
    coeff_m = 0.5 * (coeff[j] + coeff[j - 1])
    flux_p = coeff_p * (q[j + 1] - q[j]) / dz
    flux_m = coeff_m * (q[j] - q[j - 1]) / dz
    return (flux_p - flux_m) / dz
end

function radial_derivative_bdf2(current, prev1, prev2::Union{Float64,Nothing}, dr)
    if prev2 === nothing
        return (current - prev1) / dr
    end
    return (3.0 * current - 4.0 * prev1 + prev2) / (2.0 * dr)
end

function residual_station!(
    res::Vector{Float64},
    x::Vector{Float64},
    previous::Union{Station,Nothing},
    previous2::Union{Station,Nothing},
    r::Float64,
    rprev::Float64,
    z::Vector{Float64},
    p::Params,
)
    Nz = p.Nz
    dz = z[2] - z[1]
    dr = r - rprev
    U, V, W, T = unpack(x, Nz)

    fill!(res, 0.0)

    mu = [suth_mu(T[j], p) for j in 1:Nz]
    kappa = [suth_k(T[j], p) for j in 1:Nz]
    rh = [suth_rho(T[j]) for j in 1:Nz]

    # U momentum with wall/far-field boundary conditions.
    res[1] = U[1]
    for j in 2:Nz-1
        Uz = derivative_center(U, j, dz)
        if previous === nothing
            Ur = 0.0
        else
            prev2_U = previous2 === nothing ? nothing : previous2.U[j]
            Ur = radial_derivative_bdf2(U[j], previous.U[j], prev2_U, dr)
        end
        lhs = diffusive_flux_second(U, mu, j, dz)
        rhs = rh[j] * (r * U[j] * Ur + U[j]^2 + W[j] * Uz - (V[j] + 1.0)^2)
        res[j] = lhs - rhs
    end
    res[Nz] = U[Nz]

    # V momentum with wall/far-field boundary conditions.
    offV = Nz
    res[offV + 1] = V[1]
    for j in 2:Nz-1
        Vz = derivative_center(V, j, dz)
        if previous === nothing
            Vr = 0.0
        else
            prev2_V = previous2 === nothing ? nothing : previous2.V[j]
            Vr = radial_derivative_bdf2(V[j], previous.V[j], prev2_V, dr)
        end
        lhs = diffusive_flux_second(V, mu, j, dz)
        rhs = rh[j] * (r * U[j] * Vr + W[j] * Vz + 2.0 * U[j] * (V[j] + 1.0))
        res[offV + j] = lhs - rhs
    end
    res[offV + Nz] = V[Nz] + 1.0

    # Continuity determines W. For the first station, use local-similar
    # radial divergence d(r^2 rho U)/(r dr) -> 2 rho U.
    offW = 2Nz
    res[offW + 1] = W[1]
    for j in 2:Nz
        dz_flux = (rh[j] * W[j] - rh[j - 1] * W[j - 1]) / dz
        radial_div = if previous === nothing
            2.0 * rh[j] * U[j]
        else
            prev_rho = suth_rho(previous.T[j])
            flux_now = r^2 * rh[j] * U[j]
            flux_prev = rprev^2 * prev_rho * previous.U[j]
            if previous2 === nothing
                (flux_now - flux_prev) / (r * dr)
            else
                prev2_rho = suth_rho(previous2.T[j])
                flux_prev2 = previous2.r^2 * prev2_rho * previous2.U[j]
                (3.0 * flux_now - 4.0 * flux_prev + flux_prev2) / (2.0 * r * dr)
            end
        end
        res[offW + j] = dz_flux + radial_div
    end

    # Energy equation. Conductivity is k=mu/sigma, kept as variable. The
    # viscous source uses local Mach Mr(r)=Minf*r.
    offT = 3Nz
    res[offT + 1] = T[1] - p.Tw
    Mr = p.Minf * r
    for j in 2:Nz-1
        Tz = derivative_center(T, j, dz)
        Uz = derivative_center(U, j, dz)
        Vz = derivative_center(V, j, dz)
        if previous === nothing
            Tr = 0.0
        else
            prev2_T = previous2 === nothing ? nothing : previous2.T[j]
            Tr = radial_derivative_bdf2(T[j], previous.T[j], prev2_T, dr)
        end
        lhs = diffusive_flux_second(T, kappa, j, dz)
        convection = rh[j] * (r * U[j] * Tr + W[j] * Tz)
        dissipation = (p.gamma - 1.0) * Mr^2 * mu[j] * (Uz^2 + Vz^2)
        res[offT + j] = lhs - convection + dissipation
    end
    res[offT + Nz] = T[Nz] - 1.0

    # Penalize unphysical trial states during Newton line-search.
    for j in 1:Nz
        if T[j] <= 0.0 || !isfinite(T[j])
            res[offT + j] += 1.0e6 * (1.0 + abs(T[j]))
        end
    end

    return res
end

function residual_norm(res)
    return norm(res, Inf)
end

function finite_difference_jacobian(resfun, x, fx, p::Params)
    n = length(x)
    J = Matrix{Float64}(undef, n, n)
    xt = copy(x)
    ft = similar(fx)
    for k in 1:n
        old = xt[k]
        h = p.fd_eps * max(1.0, abs(old))
        xt[k] = old + h
        resfun(ft, xt)
        @views J[:, k] .= (ft .- fx) ./ h
        xt[k] = old
    end
    return J
end

function newton_solve(x0, resfun, p::Params; label="")
    x = copy(x0)
    fx = zeros(length(x))
    resfun(fx, x)
    nrm = residual_norm(fx)
    if p.verbose
        @printf("  %-18s iter=%2d  ||R||_inf=%.4e\n", label, 0, nrm)
    end
    success = nrm < p.newton_tol
    iter_done = 0

    for iter in 1:p.newton_maxiter
        success && break
        J = finite_difference_jacobian(resfun, x, fx, p)
        dx = -(J \ fx)

        accepted = false
        alpha = 1.0
        best_x = x
        best_fx = fx
        best_nrm = nrm
        for _ in 1:12
            trial = x .+ alpha .* dx
            # Keep temperature positive during trial.
            _, _, _, Ttrial = unpack(trial, p.Nz)
            if minimum(Ttrial) <= 0.0
                alpha *= 0.5
                continue
            end
            ftrial = similar(fx)
            resfun(ftrial, trial)
            trial_nrm = residual_norm(ftrial)
            if isfinite(trial_nrm) && trial_nrm < best_nrm
                best_x = trial
                best_fx = ftrial
                best_nrm = trial_nrm
                accepted = true
                break
            end
            alpha *= 0.5
        end

        x = copy(best_x)
        fx = copy(best_fx)
        nrm = best_nrm
        iter_done = iter
        if p.verbose
            @printf("  %-18s iter=%2d  ||R||_inf=%.4e  alpha=%.3g\n", label, iter, nrm, alpha)
        end
        success = accepted && nrm < p.newton_tol
        if !accepted
            break
        end
    end

    return x, success, nrm, iter_done
end

function station_from_x(r, x, p::Params)
    Uv, Vv, Wv, Tv = unpack(x, p.Nz)
    return Station(r, collect(Uv), collect(Vv), collect(Wv), collect(Tv))
end

function solve_marching(p::Params)
    z = collect(range(0.0, p.zmax, length=p.Nz))
    r_values = collect(range(p.r0, p.rmax, length=p.Nr))

    stations = Station[]

    # First station: local-similar closure for radial divergence and no radial
    # derivative terms. This supplies a physically consistent initial profile
    # for the downstream non-similar march.
    x0 = initial_guess(z, p)
    res0 = (res, x) -> residual_station!(res, x, nothing, nothing, r_values[1], r_values[1], z, p)
    x, success, nrm, _ = newton_solve(x0, res0, p; label=@sprintf("r=%.4g", r_values[1]))
    if !success
        error("Initial station failed at r=$(r_values[1]) with residual $nrm")
    end
    push!(stations, station_from_x(r_values[1], x, p))

    for i in 2:length(r_values)
        r = r_values[i]
        rprev = r_values[i - 1]
        prev = stations[end]
        prev2 = i >= 3 ? stations[end - 1] : nothing
        x_guess = pack(prev.U, prev.V, prev.W, prev.T)
        resi = (res, xx) -> residual_station!(res, xx, prev, prev2, r, rprev, z, p)
        xi, ok, rn, _ = newton_solve(x_guess, resi, p; label=@sprintf("r=%.4g", r))
        if !ok
            error("Marching failed at r=$r with residual $rn. Try smaller radial step or larger zmax.")
        end
        push!(stations, station_from_x(r, xi, p))
    end

    return r_values, z, stations
end

function write_outputs(r_values, z, stations, p::Params)
    mkpath(p.out_dir)
    profile_path = joinpath(p.out_dir, "sutherland_marching_profiles.csv")
    target_profile_path = joinpath(p.out_dir, "sutherland_profile_at_R.csv")
    summary_path = joinpath(p.out_dir, "sutherland_marching_summary.csv")
    meta_path = joinpath(p.out_dir, "metadata.txt")

    open(profile_path, "w") do io
        println(io, "r,z,U,V,W,T,rho,mu,kappa,Mr")
        for st in stations
            Mr = p.Minf * st.r
            for j in eachindex(z)
                Tj = st.T[j]
                muj = suth_mu(Tj, p)
                @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                    st.r, z[j], st.U[j], st.V[j], st.W[j], Tj, suth_rho(Tj), muj, muj / p.sigma, Mr)
            end
        end
    end

    open(target_profile_path, "w") do io
        println(io, "r,z,U,V,W,T,rho,mu,kappa,Mr,u_physical,v_rotating_physical,v_inertial_physical")
        st = stations[end]
        Mr = p.Minf * st.r
        for j in eachindex(z)
            Tj = st.T[j]
            muj = suth_mu(Tj, p)
            u_phys = st.r * st.U[j]
            v_rot_phys = st.r * st.V[j]
            v_inertial_phys = st.r * (st.V[j] + 1.0)
            @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                st.r, z[j], st.U[j], st.V[j], st.W[j], Tj, suth_rho(Tj), muj, muj / p.sigma, Mr,
                u_phys, v_rot_phys, v_inertial_phys)
        end
    end

    open(summary_path, "w") do io
        println(io, "r,Mr,Umax,z_Umax,Vmin,z_Vmin,Wmin,z_Wmin,Tmin,Tmax")
        for st in stations
            iu = argmax(st.U)
            iv = argmin(st.V)
            iw = argmin(st.W)
            @printf(io, "%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n",
                st.r, p.Minf * st.r,
                st.U[iu], z[iu],
                st.V[iv], z[iv],
                st.W[iw], z[iw],
                minimum(st.T), maximum(st.T))
        end
    end

    open(meta_path, "w") do io
        println(io, "Generated: $(Dates.now())")
        println(io, "Model: non-similar compressible rotating-disk boundary-layer marching with Sutherland viscosity")
        println(io, "Tw=$(p.Tw)")
        println(io, "sigma=$(p.sigma)")
        println(io, "gamma=$(p.gamma)")
        println(io, "Minf=$(p.Minf)")
        println(io, "Tinf_dim=$(p.Tinf_dim)")
        println(io, "S_dim=$(p.S_dim)")
        println(io, "r0=$(p.r0)")
        println(io, "rmax=$(p.rmax)")
        println(io, "final_R=$(stations[end].r)")
        println(io, "final_Mr=$(p.Minf * stations[end].r)")
        println(io, "Nr=$(p.Nr)")
        println(io, "zmax=$(p.zmax)")
        println(io, "Nz=$(p.Nz)")
        println(io, "radial_discretization=BDF2 after the second station; first station local-similar; second station backward Euler")
        println(io, "velocity_form=u=r*U, v_rotating=r*V, v_inertial=r*(V+1), w=W")
        println(io, "newton_tol=$(p.newton_tol)")
    end

    return profile_path, target_profile_path, summary_path, meta_path
end

function parse_bool(s)
    return lowercase(s) in ("1", "true", "yes", "y")
end

function usage()
    return """
Sutherland marching solver usage:

  julia SutherlandMarching.jl --Tw=1.5 --Mr=0.3 --R=20 --Nr=41 --Nz=81 --zmax=20 --out_dir=sutherland_case

Main arguments:

  --Tw=VALUE          Wall temperature, nondimensionalized by T_inf.
  --Mr=VALUE          Target local Mach number at R. If this is given, Minf is set to Mr/R.
  --R=VALUE           Target radius. This is an alias for the final marching radius rmax.
  --Minf=VALUE        Far-field/radial Mach coefficient. Local Mach is Mr(r)=Minf*r.
  --r0=VALUE          Initial marching radius.
  --rmax=VALUE        Final marching radius. Used when --R is not supplied.
  --Nr=VALUE          Number of radial stations from r0 to R/rmax.
  --Nz=VALUE          Number of wall-normal grid points.
  --zmax=VALUE        Wall-normal domain height.
  --out_dir=PATH      Output directory.

Outputs:

  sutherland_profile_at_R.csv       Profile at the requested final radius R.
  sutherland_marching_profiles.csv  Profiles at all radial stations.
  sutherland_marching_summary.csv   Peak/minimum summary at all radial stations.
  metadata.txt                      Run parameters and velocity definitions.
"""
end

function parse_args(args)
    if any(a -> a in ("--help", "-h"), args)
        println(usage())
        exit(0)
    end

    vals = Dict{String,String}()
    for a in args
        if startswith(a, "--") && occursin("=", a)
            key, val = split(a[3:end], "=", limit=2)
            vals[key] = val
        end
    end

    tw = parse(Float64, get(vals, "Tw", "1.5"))
    sigma = parse(Float64, get(vals, "sigma", "0.72"))
    gamma = parse(Float64, get(vals, "gamma", "1.4"))
    t_inf_dim = parse(Float64, get(vals, "Tinf_dim", "273.0"))
    s_dim = parse(Float64, get(vals, "S_dim", "114.0"))
    r0 = parse(Float64, get(vals, "r0", "1.0"))
    rmax = parse(Float64, get(vals, "R", get(vals, "rmax", "100.0")))
    nr = parse(Int, get(vals, "Nr", "41"))
    zmax = parse(Float64, get(vals, "zmax", "20.0"))
    nz = parse(Int, get(vals, "Nz", "81"))
    minf = parse(Float64, get(vals, "Minf", "0.003"))

    if haskey(vals, "Mr")
        target_mr = parse(Float64, vals["Mr"])
        rmax <= 0.0 && error("--R/--rmax must be positive when --Mr is supplied")
        minf = target_mr / rmax
    end

    rmax < r0 && error("The final radius R/rmax must be greater than or equal to r0")
    if isapprox(rmax, r0; atol=1.0e-12, rtol=1.0e-12)
        nr = 1
    elseif nr < 2
        error("Nr must be at least 2 when R/rmax is greater than r0")
    end

    return Params(
        Tw = tw,
        sigma = sigma,
        gamma = gamma,
        Minf = minf,
        Tinf_dim = t_inf_dim,
        S_dim = s_dim,
        r0 = r0,
        rmax = rmax,
        Nr = nr,
        zmax = zmax,
        Nz = nz,
        newton_tol = parse(Float64, get(vals, "newton_tol", "1e-8")),
        newton_maxiter = parse(Int, get(vals, "newton_maxiter", "14")),
        fd_eps = parse(Float64, get(vals, "fd_eps", "1e-6")),
        out_dir = get(vals, "out_dir", "sutherland_marching_data"),
        verbose = parse_bool(get(vals, "verbose", "true")),
    )
end

function main()
    p = parse_args(ARGS)
    @printf("Sutherland marching solver: Tw=%.4g, Minf=%.4g, R=%.4g, Mr(R)=%.4g, r=[%.4g, %.4g], Nr=%d, Nz=%d\n",
        p.Tw, p.Minf, p.rmax, p.Minf * p.rmax, p.r0, p.rmax, p.Nr, p.Nz)
    r_values, z, stations = solve_marching(p)
    profile_path, target_profile_path, summary_path, meta_path = write_outputs(r_values, z, stations, p)
    println("Done.")
    println("Profile at R: ", target_profile_path)
    println("All profiles: ", profile_path)
    println("Summary:      ", summary_path)
    println("Metadata:     ", meta_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module SutherlandMarching

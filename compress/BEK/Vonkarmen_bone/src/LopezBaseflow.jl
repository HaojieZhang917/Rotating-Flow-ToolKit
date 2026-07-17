module LopezBaseflow

using DelimitedFiles
using Printf
using PyCall

export solve_lopez_baseflow, get_baseflow, write_profile, clear_cache!

const DEFAULT_PR = 0.72
const PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const PYTHON_MODULE = Ref{Union{Nothing, PyObject}}(nothing)
const PROFILE_CACHE = Dict{Tuple, Any}()

function python_module()
    if PYTHON_MODULE[] === nothing
        sys = pyimport("sys")
        pushfirst!(PyVector(sys."path"), PROJECT_ROOT)
        PYTHON_MODULE[] = pyimport("compare_lopez_boussinesq")
    end
    return PYTHON_MODULE[]::PyObject
end

"""Clear profiles cached by `solve_lopez_baseflow`."""
function clear_cache!()
    empty!(PROFILE_CACHE)
    return nothing
end

function scipy_field(solution, name::AbstractString)
    if solution isa AbstractDict
        return solution[name]
    end
    return getproperty(solution, name)
end

function profile_metrics(Tw, z, H, F, G, T, Fp, Gp, Tp, solution, step)
    peak_index = argmax(F)
    boundary_residuals = [
        H[1], F[1], G[1], T[1] - Tw,
        F[end], G[end] - 1.0, T[end] - 1.0,
    ]
    return (
        Hinf=H[end],
        Fmax=F[peak_index],
        z_Fmax=z[peak_index],
        Fp0=Fp[1],
        Gp0=Gp[1],
        Tp0=Tp[1],
        rho_wall_linear=2.0 - Tw,
        boundary_residual_norm=maximum(abs, boundary_residuals),
        solver_nodes=length(scipy_field(solution, "x")),
        continuation_steps=ceil(Int, abs(Tw - 1.0) / step),
        success=Bool(scipy_field(solution, "success")),
        message=String(scipy_field(solution, "message")),
        backend="scipy.integrate.solve_bvp via PyCall",
    )
end

"""
    solve_lopez_baseflow(Tw; kwargs...)

Call the validated SciPy Lopez base-flow solver through PyCall and return a
Julia named tuple containing `z,H,F,G,T,Fp,Gp,Tp,metrics,raw_solution`.

Keywords:

- `Pr=0.72`: currently fixed by the Python model.
- `zmax=40.0`: finite far-field boundary.
- `continuation_step=0.0025`: maximum wall-temperature continuation step.
- `tol=1e-8`: SciPy BVP tolerance.
- `nout=2001`: number of uniformly sampled output points.
- `cache=true`: reuse an identical solution in later notebook cells.
- `save=false`: do not write a file unless explicitly requested.
- `output_path=nothing`: optional CSV path when `save=true`.
- `verbose=false`: suppress Python continuation progress.

`G` is the Lopez azimuthal deficit with `G(0)=0`, `G(infinity)=1`.
"""
function solve_lopez_baseflow(
    Tw::Real;
    Pr::Real=DEFAULT_PR,
    zmax::Real=40.0,
    continuation_step::Real=0.0025,
    tol::Real=1.0e-8,
    nout::Integer=2001,
    cache::Bool=true,
    save::Bool=false,
    output_path::Union{Nothing, AbstractString}=nothing,
    verbose::Bool=false,
)
    Tw >= 1.0 || throw(ArgumentError(
        "the existing Python continuation currently supports Tw >= 1",
    ))
    Tw < 2.0 || throw(ArgumentError(
        "Tw must remain below 2 so the linear density approximation is positive",
    ))
    isapprox(Pr, DEFAULT_PR; atol=10eps(Float64), rtol=0.0) || throw(ArgumentError(
        "the existing Python solver is fixed at Pr=0.72",
    ))
    zmax > 0 || throw(ArgumentError("zmax must be positive"))
    continuation_step > 0 || throw(ArgumentError("continuation_step must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    nout >= 2 || throw(ArgumentError("nout must be at least 2"))

    Tw64 = Float64(Tw)
    zmax64 = Float64(zmax)
    step64 = Float64(continuation_step)
    tol64 = Float64(tol)
    key = (Tw64, Float64(Pr), zmax64, step64, tol64, Int(nout))

    if cache && haskey(PROFILE_CACHE, key)
        profile = PROFILE_CACHE[key]
    else
        solver = python_module()
        _, solutions = solver.continue_to_targets(
            "lopez", [Tw64];
            zmax=zmax64, step=step64, tol=tol64, verbose=verbose,
        )
        solution = solutions[Tw64]
        Bool(scipy_field(solution, "success")) || error(
            "SciPy Lopez solve failed at Tw=$Tw64: " *
            String(scipy_field(solution, "message")),
        )

        z = collect(range(0.0, zmax64; length=Int(nout)))
        state = Matrix{Float64}(scipy_field(solution, "sol")(z))

        H = vec(state[1, :])
        Fp = vec(state[2, :])
        F = vec(state[3, :])
        Gp = vec(state[4, :])
        G = vec(state[5, :])
        Tp = vec(state[6, :])
        T = vec(state[7, :])
        metrics = profile_metrics(
            Tw64, z, H, F, G, T, Fp, Gp, Tp, solution, step64,
        )
        profile = (
            Tw=Tw64, Pr=Float64(Pr), zmax=zmax64,
            z=z, H=H, F=F, G=G, T=T, Fp=Fp, Gp=Gp, Tp=Tp,
            metrics=metrics, raw_solution=solution,
        )
        cache && (PROFILE_CACHE[key] = profile)
    end

    if save
        path = output_path === nothing ? default_output_path(Tw64) : String(output_path)
        write_profile(path, profile)
    end
    return profile
end

"""
    get_baseflow(Tw; kwargs...)

Compatibility interface matching the current Python `bone.get_baseflow`
layout. Returns `z,H,F,G,T,Fp,Gp,Tp,info` as Julia values.
"""
function get_baseflow(Tw::Real; kwargs...)
    profile = solve_lopez_baseflow(Tw; kwargs...)
    return (
        profile.z, profile.H, profile.F, profile.G, profile.T,
        profile.Fp, profile.Gp, profile.Tp, profile.metrics,
    )
end

function default_output_path(Tw)
    tag = replace(@sprintf("%.6g", Tw), "." => "p", "-" => "m")
    return "lopez_baseflow_Tw$(tag).csv"
end

"""Write a profile using columns compatible with the stability routines."""
function write_profile(path::AbstractString, profile)
    parent = dirname(path)
    parent == "." || mkpath(parent)
    data = hcat(
        profile.z, profile.H, profile.F, profile.G, profile.T,
        profile.Fp, profile.Gp, profile.Tp,
    )
    open(path, "w") do stream
        println(stream, "z,H,F,G,T,Fp,Gp,Tp")
        writedlm(stream, data, ',')
    end
    return abspath(path)
end

end # module LopezBaseflow

#!/usr/bin/env julia

"""Adjoint projection analysis of the traditional Boussinesq fold.

The calculation is performed at the first, isothermal-connected fold for
Re_h=1000 and Pr=0.72.  The fixed-Tw base-flow Jacobian is decomposed into
right/left near-null vectors.  We report (i) projections of the residual
terms, which describe the force balance at the fold, and (ii) projections of
the quadratic terms along the right null vector, which diagnose the
nonlinearity that creates the saddle-node curvature.

This is an adjoint analysis of the steady basic-flow Jacobian.  It is not the
adjoint of the perturbation LST operator.
"""

using CompressRotatingCavity
using LinearAlgebra

const TB = CompressRotatingCavity.TraditionalBoussinesq
const ROOT = normpath(joinpath(@__DIR__, ".."))
const REFERENCE = joinpath(ROOT, "data", "reference", "baseflow_Res1000.npz")
const OUTPUT = joinpath(ROOT, "benchmark_results")

function residual_terms(profile)
    n = length(profile.z)
    m = 4n + 1
    D = profile.D
    D2 = D * D
    fields = profile.fields
    H = vec(fields[1, :])
    F = vec(fields[2, :])
    G = vec(fields[3, :])
    T = vec(fields[4, :])
    cfg = profile.config
    re = cfg.re_h
    sqrt_re = sqrt(re)
    pr = cfg.pr
    DH, DF, DG, DT = D * H, D * F, D * G, D * T

    names = [
        "H_diffusion", "H_F_coupling",
        "F_diffusion", "F_axial_advection", "F_radial_inertia",
        "F_rotation", "F_thermal", "F_pressure",
        "G_diffusion", "G_axial_advection", "G_FG_coupling",
        "T_diffusion", "T_axial_advection",
    ]
    vectors = Dict(name => zeros(Float64, m) for name in names)

    rH = 1:n
    rF = n+1:2n
    rG = 2n+1:3n
    rT = 3n+1:4n
    # Boundary rows are deliberately left zero.  The pressure and wall
    # temperature sensitivities are added separately below.
    iH = 2:n
    iF = (n + 2):(2n - 1)
    iG = (2n + 2):(3n - 1)
    iT = (3n + 2):(4n - 1)

    vectors["H_diffusion"][iH] .= DH[2:n]
    vectors["H_F_coupling"][iH] .= (2sqrt_re .* F)[2:n]

    vectors["F_diffusion"][iF] .= (D2 * F)[2:n-1]
    vectors["F_axial_advection"][iF] .= (-sqrt_re .* H .* DF)[2:n-1]
    vectors["F_radial_inertia"][iF] .= (-re .* F.^2)[2:n-1]
    vectors["F_rotation"][iF] .= (re .* G.^2)[2:n-1]
    vectors["F_thermal"][iF] .= (-re .* (T .- 1.0))[2:n-1]
    vectors["F_pressure"][iF] .= -re .* profile.pressure

    vectors["G_diffusion"][iG] .= (D2 * G)[2:n-1]
    vectors["G_axial_advection"][iG] .= (-sqrt_re .* H .* DG)[2:n-1]
    vectors["G_FG_coupling"][iG] .= (-2re .* F .* G)[2:n-1]

    vectors["T_diffusion"][iT] .= (D2 * T)[2:n-1]
    vectors["T_axial_advection"][iT] .=
        (-pr * sqrt_re .* H .* DT)[2:n-1]
    return vectors
end

function quadratic_terms(profile, phi)
    n = length(profile.z)
    m = 4n + 1
    D = profile.D
    cfg = profile.config
    re = cfg.re_h
    sqrt_re = sqrt(re)
    pr = cfg.pr
    phiH = view(phi, 1:n)
    phiF = view(phi, n+1:2n)
    phiG = view(phi, 2n+1:3n)
    phiT = view(phi, 3n+1:4n)
    dphiF = D * phiF
    dphiG = D * phiG
    dphiT = D * phiT

    names = [
        "F_HDF_quadratic", "F_F2_quadratic", "F_G2_quadratic",
        "G_HDG_quadratic", "G_FG_quadratic", "T_HDT_quadratic",
    ]
    vectors = Dict(name => zeros(Float64, m) for name in names)
    iF = (n + 2):(2n - 1)
    iG = (2n + 2):(3n - 1)
    iT = (3n + 2):(4n - 1)

    vectors["F_HDF_quadratic"][iF] .=
        (-2sqrt_re .* phiH .* dphiF)[2:n-1]
    vectors["F_F2_quadratic"][iF] .= (-2re .* phiF.^2)[2:n-1]
    vectors["F_G2_quadratic"][iF] .= (2re .* phiG.^2)[2:n-1]
    vectors["G_HDG_quadratic"][iG] .=
        (-2sqrt_re .* phiH .* dphiG)[2:n-1]
    vectors["G_FG_quadratic"][iG] .= (-4re .* phiF .* phiG)[2:n-1]
    vectors["T_HDT_quadratic"][iT] .=
        (-2pr * sqrt_re .* phiH .* dphiT)[2:n-1]
    return vectors
end

function parameter_vectors(profile)
    n = length(profile.z)
    m = 4n + 1
    re = profile.config.re_h
    dPi = zeros(Float64, m)
    dTw = zeros(Float64, m)
    iF = (n + 2):(2n - 1)
    dPi[iF] .= -re
    # The first temperature row is T(0)-Tw=0 in the fixed-temperature system.
    dTw[3n + 1] = -1.0
    return dPi, dTw
end

function projection(adjoint, vector)
    return real(dot(adjoint, vector))
end

function write_vector_table(path, z, phi, psi)
    n = length(z)
    open(path, "w") do io
        println(io, "z,phi_H,phi_F,phi_G,phi_T,psi_H,psi_F,psi_G,psi_T")
        for k in 1:n
            println(io, join((z[k], phi[k], phi[n+k], phi[2n+k], phi[3n+k],
                              psi[k], psi[n+k], psi[2n+k], psi[3n+k]), ','))
        end
    end
end

function main()
    degree = isempty(ARGS) ? 80 : parse(Int, ARGS[1])
    mkpath(OUTPUT)
    cfg = ThermalConfig(re_h=1000.0, pr=0.72, degree=degree,
                        tolerance=2e-8)
    iso = isothermal_profile(cfg, REFERENCE)
    fold, _ = trace_principal_fold(iso; coarse_step=0.003,
                                   fine_step=0.0005,
                                   pressure_tolerance=2e-9)

    state = vcat(TB.pack(fold.fields), fold.pressure)
    residual, J = TB.system(state, cfg, :temperature, fold.tw, fold.D)
    decomposition = svd(Matrix(J))
    sigma_min = decomposition.S[end]
    phi = copy(decomposition.Vt[end, :])
    psi = copy(decomposition.U[:, end])
    phi ./= max(norm(phi[1:4length(fold.z)]), eps())
    psi ./= max(norm(psi), eps())

    # Fix the arbitrary signs for reproducible tables.  The sign has no effect
    # on normalized term magnitudes, but makes profiles comparable across N.
    pivot = argmax(abs.(phi[1:4length(fold.z)]))
    phi[pivot] < 0 && (phi .*= -1)
    dot(psi, phi) < 0 && (psi .*= -1)

    terms = residual_terms(fold)
    quadratic = quadratic_terms(fold, phi)
    dPi, dTw = parameter_vectors(fold)
    pressure_projection = projection(psi, dPi)
    temperature_projection = projection(psi, dTw)
    residual_norm = norm(residual, Inf)
    h_fd = 1e-3
    residual_plus = first(TB.system(state .+ h_fd .* phi, cfg,
                                    :temperature, fold.tw, fold.D))
    residual_minus = first(TB.system(state .- h_fd .* phi, cfg,
                                     :temperature, fold.tw, fold.D))
    quadratic_fd = projection(
        psi,
        (residual_plus .+ residual_minus .- 2.0 .* residual) ./ h_fd^2,
    )

    term_values = [("residual", name, projection(psi, vector))
                   for (name, vector) in sort(collect(terms); by=first)]
    quadratic_values = [("quadratic", name, projection(psi, vector))
                        for (name, vector) in sort(collect(quadratic); by=first)]
    all_values = vcat(term_values, quadratic_values)
    scales = Dict("residual" => sum(abs(row[3]) for row in term_values),
                  "quadratic" => sum(abs(row[3]) for row in quadratic_values))
    residual_projection_sum = sum(row[3] for row in term_values)
    quadratic_projection_sum = sum(row[3] for row in quadratic_values)

    suffix = degree == 80 ? "" : "_N$(degree)"
    term_path = joinpath(OUTPUT, "fold_projection_terms$(suffix).csv")
    open(term_path, "w") do io
        println(io, "kind,name,projection,abs_fraction")
        for (kind, name, value) in all_values
            denom = max(scales[kind], eps())
            println(io, join((kind, name, value, abs(value) / denom), ','))
        end
    end

    vector_path = joinpath(OUTPUT, "fold_projection_modes$(suffix).csv")
    write_vector_table(vector_path, fold.z, phi, psi)

    summary_path = joinpath(OUTPUT, "fold_projection_summary$(suffix).txt")
    open(summary_path, "w") do io
        println(io, "Adjoint projection at the traditional Boussinesq fold")
        println(io, "Re_h=$(cfg.re_h), Pr=$(cfg.pr), degree=$(degree)")
        println(io, "Tw_fold=$(fold.tw), Pi_fold=$(fold.pressure)")
        println(io, "baseflow_residual_inf=$(residual_norm)")
        println(io, "smallest_J_singular_value=$(sigma_min)")
        println(io, "dot(adjoint,right_mode)=$(dot(psi, phi))")
        println(io, "projection_dR_dPi_augmented_column=$(pressure_projection)")
        println(io, "projection_dR_dTw=$(temperature_projection)")
        println(io, "residual_projection_sum=$(residual_projection_sum)")
        println(io, "quadratic_projection_sum=$(quadratic_projection_sum)")
        println(io, "quadratic_projection_finite_difference=$(quadratic_fd)")
        println(io, "")
        println(io, "Residual-balance projections")
        println(io, "name,projection,abs_fraction")
        for (_, name, value) in term_values
            println(io, join((name, value,
                              abs(value) / max(scales["residual"], eps())), ','))
        end
        println(io, "")
        println(io, "Quadratic fold-curvature projections")
        println(io, "name,projection,abs_fraction")
        for (_, name, value) in quadratic_values
            println(io, join((name, value,
                              abs(value) / max(scales["quadratic"], eps())), ','))
        end
        println(io, "")
        println(io, "The Pi projection is near zero because Pi is an unknown column")
        println(io, "of the augmented fixed-Tw Jacobian at the fold.")
        println(io, "The quadratic values are projections of R_yy[phi,phi].")
        println(io, "They diagnose the nonlinear curvature coefficient of the fold;")
        println(io, "they are distinct from perturbation-LST adjoint projections.")
    end

    println("Tw_fold=$(fold.tw), Pi_fold=$(fold.pressure)")
    println("smallest J singular value=$(sigma_min)")
    println("projection dR/dPi=$(pressure_projection)")
    println("projection dR/dTw=$(temperature_projection)")
    println("quadratic projection (analytic)=$(quadratic_projection_sum)")
    println("quadratic projection (finite difference)=$(quadratic_fd)")
    println("Wrote ", term_path)
    println("Wrote ", vector_path)
    println("Wrote ", summary_path)
end

main()

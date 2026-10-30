#!/usr/bin/env julia

"""Local steady axisymmetric spatial spectrum of the single rotating disk.

Run `julia --project=. scripts/spectrum_single_disk_radial.jl`.
The pencil is J + kappa*B for perturbations exp(kappa*s), s=log(r/r0).
This is not a temporal or three-dimensional stability calculation.
Only Julia standard libraries and the directly included radial solver are used.
"""

using LinearAlgebra

include(joinpath(@__DIR__, "..", "src", "SingleDiskRadial.jl"))
using .SingleDiskRadial

const SPECTRUM_OUTPUT = normpath(joinpath(@__DIR__, "..", "benchmark_results", "single_disk_radial"))
const SPECTRAL_WINDOW = 50.0
const EIGEN_RESIDUAL_TOLERANCE = 1e-10
const GRID_RELATIVE_TOLERANCE = 1e-7
const DOMAIN_RELATIVE_TOLERANCE = 1e-5
const MASS_COEFFICIENT_TOLERANCE = 1e-10
const MODE_CORRELATION_TOLERANCE = 0.995

function csv_cell(value)
    text = string(value)
    return occursin(r"[,\"\n\r]", text) ? "\"" * replace(text, "\"" => "\"\"") * "\"" : text
end

"""Remove both the prescribed perturbation DOFs and their boundary rows.

Continuity contributes +kappa*f, and the two momentum equations contribute
-kappa*F0*f and -kappa*F0*g.  H has no independent radial derivative.
The remaining singular mass matrix has diagnostic/infinite eigenvalues.
"""
function spatial_pencil(reference)
    n = length(reference.eta)
    _, J = radial_system(reference.state, reference)
    B = zeros(3n, 3n)
    B[1:n, (n+1):(2n)] .= Matrix{Float64}(I, n, n)
    B[(n+1):(2n), (n+1):(2n)] .= -Diagonal(reference.F)
    B[(2n+1):(3n), (2n+1):(3n)] .= -Diagonal(reference.F)
    prescribed = [1, n+1, 2n, 2n+1, 3n]
    B[prescribed, :] .= 0
    free = setdiff(collect(1:3n), prescribed)

    # Keep Fs=Gs=0 at the base while differentiating their dependence on F,G.
    # With histories fixed in this derivative, J(a=1)-J(a=0) is exactly B.
    _, J1 = radial_system(reference.state, reference; a=1.0,
                          history_f=-reference.F, history_g=-reference.G)
    mass_error = norm(J1-J-B, Inf) / max(norm(B, Inf), eps(Float64))
    return (; J=J[free, free], B=B[free, free], free, mass_error)
end

function solve_spectrum(; degree, eta_max)
    config = RadialConfig(; degree, eta_max)
    reference = reference_profile(config)
    pencil = spatial_pencil(reference)
    factorization = eigen(pencil.J, -pencil.B)
    normJ, normB = opnorm(pencil.J, Inf), opnorm(pencil.B, Inf)
    modes = NamedTuple[]
    for i in eachindex(factorization.values)
        kappa = ComplexF64(factorization.values[i])
        finite = isfinite(real(kappa)) && isfinite(imag(kappa))
        vector = ComplexF64.(factorization.vectors[:, i])
        denominator = (normJ+abs(kappa)*normB)*norm(vector, Inf)
        residual = finite && isfinite(denominator) && denominator > 0 ?
            norm(pencil.J*vector+kappa*(pencil.B*vector), Inf)/denominator : Inf
        in_window = finite && abs(kappa) <= SPECTRAL_WINDOW
        admissible = in_window && residual <= EIGEN_RESIDUAL_TOLERANCE
        full = zeros(ComplexF64, length(reference.state))
        full[pencil.free] .= vector
        push!(modes, (; raw_index=i, kappa, finite, residual, in_window, admissible,
                      full, constrained_dof_error=maximum(abs.(full[setdiff(1:length(full), pencil.free)]))))
    end
    sort!(modes; by=mode -> mode.finite ? (-real(mode.kappa), abs(imag(mode.kappa))) : (Inf, Inf))
    return (; reference, pencil, modes)
end

function mode_samples(result, mode, points)
    n = length(result.reference.eta)
    parts = ComplexF64[]
    for range in ((n+1):(2n), (2n+1):(3n))
        values = mode.full[range]
        samples = interpolate_profile(result.reference.eta, real.(values), points) .+
            im .* interpolate_profile(result.reference.eta, imag.(values), points)
        append!(parts, samples)
    end
    return parts
end

function match_mode(result, target_result, target_mode)
    candidates = filter(mode -> mode.admissible, result.modes)
    isempty(candidates) && return (; difference=Inf, correlation=0.0, kappa=ComplexF64(NaN))
    index = argmin([abs(mode.kappa-target_mode.kappa) for mode in candidates])
    match = candidates[index]
    points = collect(range(0.0, min(12.0, result.reference.config.eta_max,
                                  target_result.reference.config.eta_max); length=301))
    a = mode_samples(result, match, points)
    b = mode_samples(target_result, target_mode, points)
    correlation = abs(dot(a, b)) / (norm(a)*norm(b))
    return (; difference=abs(match.kappa-target_mode.kappa), correlation, kappa=match.kappa)
end

"""An independent sign check: march a small leading eigenmode and fit its amplitude.

The inlet specifies only F,G.  The first computed H is solved by continuity.
The fitted rate must follow exp(kappa*s), not exp(-kappa*s).
"""
function marching_sign_check(result, mode)
    reference = result.reference
    n = length(reference.eta)
    abs(imag(mode.kappa)) < 1e-8 || return (; pass=false, rate=NaN, predicted=real(mode.kappa),
        status=:complex_leading_mode, amplitude=NaN)
    vector = real.(mode.full)
    fg = vector[(n+1):(3n)]
    fg ./= maximum(abs, fg)
    f, g = fg[1:n], fg[(n+1):(2n)]
    amplitude = 1e-6
    # Keep the test in the outward-flow branch, including the small far-field tail.
    for j in 2:(n-1)
        abs(f[j]) > 0 && (amplitude = min(amplitude, 0.01reference.F[j]/abs(f[j])))
    end
    inlet = (; case=:leading_eigenmode, F=reference.F+amplitude*f,
              G=reference.G+amplitude*g)
    config = RadialConfig(degree=reference.config.degree, eta_max=reference.config.eta_max,
                          ds=0.001, s_end=0.08, tolerance=2e-12)
    marched = march(reference, inlet; config)
    weights = vcat(reference.weights, reference.weights)
    perturbations = vcat(marched.F .- reference.F, marched.G .- reference.G)
    coefficients = vec(transpose(weights .* fg) * perturbations) ./ sum(weights .* fg.^2)
    useful = findall(s -> 0.02 <= s <= 0.08, marched.s)
    if marched.status != :completed || isempty(useful) || any(coefficients[useful] .<= 0)
        return (; pass=false, rate=NaN, predicted=real(mode.kappa), status=marched.status, amplitude)
    end
    ss = marched.s[useful]
    logs = log.(coefficients[useful])
    centered = ss .- sum(ss)/length(ss)
    rate = dot(centered, logs .- sum(logs)/length(logs)) / dot(centered, centered)
    passed = rate < 0 && abs(rate-real(mode.kappa)) <= 0.01max(1.0, abs(real(mode.kappa)))
    return (; pass=passed, rate, predicted=real(mode.kappa), status=marched.status, amplitude)
end

function main()
    BLAS.set_num_threads(1)
    configurations = [(48, 20.0), (64, 20.0), (80, 20.0), (80, 16.0), (80, 24.0)]
    spectra = [solve_spectrum(degree=N, eta_max=L) for (N,L) in configurations]
    target = spectra[3]
    window_modes = filter(mode -> mode.admissible, target.modes)
    comparison = NamedTuple[]
    for mode in window_modes
        grid48 = match_mode(spectra[1], target, mode)
        grid64 = match_mode(spectra[2], target, mode)
        domain16 = match_mode(spectra[4], target, mode)
        domain24 = match_mode(spectra[5], target, mode)
        scale = max(1.0, abs(mode.kappa))
        grid_converged = max(grid48.difference, grid64.difference) <= GRID_RELATIVE_TOLERANCE*scale &&
            min(grid48.correlation, grid64.correlation) >= MODE_CORRELATION_TOLERANCE
        domain_converged = domain24.difference <= DOMAIN_RELATIVE_TOLERANCE*scale &&
            min(domain16.correlation, domain24.correlation) >= MODE_CORRELATION_TOLERANCE &&
            domain24.difference <= max(domain16.difference, 1e-8)
        push!(comparison, (; mode, grid48, grid64, domain16, domain24,
                           grid_converged, domain_converged, converged=grid_converged && domain_converged))
    end
    accepted = filter(row -> row.converged, comparison)
    sort!(accepted; by=row -> -real(row.mode.kappa))
    sign_check = isempty(accepted) ? nothing : marching_sign_check(target, first(accepted).mode)
    # Subtracting large second-derivative diagonals incurs floating-point cancellation.
    mass_pass = all(result.pencil.mass_error <= MASS_COEFFICIENT_TOLERANCE for result in spectra)

    mkpath(SPECTRUM_OUTPUT)
    open(joinpath(SPECTRUM_OUTPUT, "spectrum.csv"), "w") do io
        println(io, "degree,eta_max,raw_index,kappa_real,kappa_imag,is_finite,scaled_eigen_residual,in_abs_kappa_le50_window,residual_accepted,constrained_dof_error,grid48_difference,grid64_difference,domain16_difference,domain24_difference,grid48_mode_correlation,grid64_mode_correlation,domain16_mode_correlation,domain24_mode_correlation,grid_converged,domain_converged,converged_physical_candidate")
        for result in spectra, mode in result.modes
            row_index = result === target ? findfirst(row -> row.mode.raw_index == mode.raw_index, comparison) : nothing
            row = row_index === nothing ? nothing : comparison[row_index]
            diagnostics = row === nothing ? fill("", 11) :
                Any[row.grid48.difference, row.grid64.difference, row.domain16.difference, row.domain24.difference,
                    row.grid48.correlation, row.grid64.correlation, row.domain16.correlation, row.domain24.correlation,
                    row.grid_converged, row.domain_converged, row.converged]
            cells = Any[result.reference.config.degree, result.reference.config.eta_max,
                        mode.raw_index, real(mode.kappa), imag(mode.kappa), mode.finite,
                        mode.residual, mode.in_window, mode.admissible, mode.constrained_dof_error]
            append!(cells, diagnostics)
            println(io, join(csv_cell.(cells), ','))
        end
    end
    open(joinpath(SPECTRUM_OUTPUT, "spectrum.txt"), "w") do io
        println(io, "Local steady axisymmetric spatial spectrum: single rotating disk main von Karman branch")
        println(io, "s=log(r/r0); perturbations exp(kappa*s)=(r/r0)^kappa. Negative real kappa means downstream decay of normalized F,G,H.")
        println(io, "Linear equations: h'+(2+kappa)f=0; f''-H0*f'-F0'*h-(2+kappa)F0*f+2G0*g=0; g''-H0*g'-G0'*h-2G0*f-(2+kappa)F0*g=0.")
        println(io, "J + kappa*B = 0; B_HF=I, B_FF=B_GG=-diag(F0), other blocks zero; generalized eigen(J,-B).")
        println(io, "State order H,F,G. Both boundary rows and DOFs Hwall,Fwall,Ftop,Gwall,Gtop removed. Htop is free.")
        println(io, "Search/report window |kappa|<=", SPECTRAL_WINDOW, "; scaled residual <=", EIGEN_RESIDUAL_TOLERANCE, ".")
        println(io, "Grid criterion: max(N48-N80,N64-N80) <= ", GRID_RELATIVE_TOLERANCE, "*max(1,|kappa|).")
        println(io, "Domain criterion: |L20-L24| <= ", DOMAIN_RELATIVE_TOLERANCE, "*max(1,|kappa|), no larger than |L16-L20| (floor 1e-8).")
        println(io, "Both comparisons also require inner-eta F,G mode correlations >=", MODE_CORRELATION_TOLERANCE, ".")
        println(io, "The singular mass matrix generates infinite diagnostic eigenvalues; finite values outside the declared window are retained in CSV but are not classified as physical.")
        println(io, "Mass/sign coefficient check: ", mass_pass ? "PASS" : "FAIL", "; relative tolerance=", MASS_COEFFICIENT_TOLERANCE, " (cancellation in J(a=1)-J(a=0)).")
        for result in spectra
            finite = count(mode -> mode.finite, result.modes)
            candidates = count(mode -> mode.admissible, result.modes)
            positive_outside = count(mode -> mode.finite && !mode.in_window && real(mode.kappa)>0, result.modes)
            println(io, "N=", result.reference.config.degree, ", L=", result.reference.config.eta_max,
                    ": total=", length(result.modes), ", finite=", finite,
                    ", residual-accepted in-window=", candidates,
                    ", unclassified positive finite outside-window=", positive_outside,
                    ", mass coefficient error=", result.pencil.mass_error)
        end
        println(io, "Some coarser discretizations can produce very large positive finite eigenvalues outside the window. These have not passed grid/domain convergence and are unresolved numerical candidates, not evidence of a physical growing mode.")
        println(io, "N80 L20 candidates and convergence diagnostics:")
        for row in comparison
            println(io, row.converged ? "CONVERGED " : "UNRESOLVED ", row.mode.kappa,
                    "; residual=", row.mode.residual,
                    "; grid delta48/64=", row.grid48.difference, "/", row.grid64.difference,
                    "; domain delta16/24=", row.domain16.difference, "/", row.domain24.difference,
                    "; grid correlations=", row.grid48.correlation, "/", row.grid64.correlation,
                    "; domain correlations=", row.domain16.correlation, "/", row.domain24.correlation)
        end
        if isempty(accepted)
            println(io, "No converged rate is claimed.")
        else
            println(io, "Leading converged candidate in the declared window: ", first(accepted).mode.kappa)
            println(io, "Number of converged in-window candidates with Re(kappa)>=0: ", count(row -> real(row.mode.kappa) >= 0, accepted))
            println(io, "Independent short nonlinear eigenmode march sign check: ", sign_check)
        end
        println(io, "kappa=-2 was not divided out: it gives h=0 and w=f+i*g obeying w''-H0*w'-2i*G0*w=0 with homogeneous endpoint values; a positive integrating-factor energy identity excludes a nonzero mode there.")
        println(io, "Scope: local steady axisymmetric spatial attraction of the normalized similarity profiles in this boundary-layer model; not temporal stability, non-axisymmetric/3D stability, global nonlinear attraction, or a proof about the full Navier-Stokes flow.")
        println(io, "The degenerating far-field radial velocity makes high-|kappa| numerical eigenvalues delicate. Convergence inside this finite window does not exclude unresolved eigenvalues elsewhere or establish a function-space spectrum theorem.")
    end
    println("Spectrum files: ", SPECTRUM_OUTPUT)
    for row in accepted
        println("Converged in-window kappa = ", row.mode.kappa)
    end
    println("Mass check: ", mass_pass, "; short-march sign check: ", sign_check)
    mass_pass || error("mass coefficient check failed")
    sign_check === nothing || sign_check.pass || error("independent marching sign check failed")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

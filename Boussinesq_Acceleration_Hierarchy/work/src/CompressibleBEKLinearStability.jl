module CompressibleBEKLinearStability

using DelimitedFiles
using LinearAlgebra
using Printf
using Random
using DifferentialEquations
using BoundaryValueDiffEq
using BSplineKit

const TOOLKIT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const BACKEND_SOURCE = joinpath(TOOLKIT_ROOT, "RotatingDiskFlow", "src", "RotatingDiskFlow.jl")
isfile(BACKEND_SOURCE) || error("RotatingDiskFlow backend not found at $BACKEND_SOURCE")
include(BACKEND_SOURCE)

const RDF = RotatingDiskFlow
const CRD = RDF.NeutralCurveRunner
const _SOURCE_PROFILE_CACHE = Dict{Tuple{Float64,Float64,Float64},Any}()

export PhysicalParameters, NumericalParameters, PreparedCase, StabilityResult
export stability_grid, prepare_case, assemble_qep, solve_mode, validate_result
export audit_expanded_operator, operator_structure_gates, directional_qep_test, write_result
export at_fixed_system_mach

"""Physical/nondimensional inputs for a local compressible BEK calculation.

`Mr` is the local differential-rotation Mach number `r*abs(DeltaOmega)/a_inf`.
For the von Karman endpoint this is the rotational Mach number used by
Turkyilmazoglu & Uygun (2006).  The primitive-equation Mach number passed to
the operator is `Ma = Mr/R`.  The published von Karman equations correspond
to `Ro = -1`.
"""
Base.@kwdef struct PhysicalParameters
    Tw::Float64 = 1.0
    Mr::Float64 = 0.1
    Ro::Float64 = -1.0
    R::Float64 = 285.36
    beta::Float64 = 0.07759
    omega::Float64 = 0.0
    gamma::Float64 = 1.4
    Pr::Float64 = 0.72
end

"""Chebyshev grid consistent with the derivative mapping on `[0, Inf)`.

The inherited rational-grid routine constructs the differentiation matrix for
an infinite algebraic map but clips every returned coordinate above `ymax` to
`ymax`.  That creates repeated far-field coordinates.  Here the same map and
the same differentiation matrices are retained, while the last coordinate is
represented explicitly by `Inf` and every finite coordinate remains unique.
"""
function stability_grid(n)
    theta = range(0.0, stop=pi, length=n.N + 1)
    eta = reshape(-cos.(theta), n.N + 1, 1)
    weights = [2; ones(n.N - 1); 2] .* (-1.0) .^ (0:n.N)
    eta_matrix = repeat(eta, 1, n.N + 1)
    delta = eta_matrix - eta_matrix'
    Deta = (weights * (1 ./ weights)') ./ (delta .+ I(n.N + 1))
    Deta -= diagm(vec(sum(Deta, dims=2)))

    if n.domain == :finite
        D = (2 / n.ymax) .* Deta
        y = (n.ymax / 2) .* (eta .+ 1)
    else
        a = 6.0
        b = 0.6
        shape = 0.5
        h = b .* eta .+ (1 - b) .* (
            eta .^ 3 .+ shape .* (1 .- eta .^ 2)
        )
        dh = b .+ (1 - b) .* (3 .* eta .^ 2 .- 2 .* shape .* eta)
        D = ((1 .- h) .^ 2 ./ (2a .* dh)) .* Deta
        denominator = 1 .- h
        y = map(eachindex(denominator)) do index
            denominator[index] <= eps(Float64) ? Inf :
                a * (1 + h[index]) / denominator[index]
        end
        y = reshape(Float64.(y), n.N + 1, 1)
    end
    return Matrix{Float64}(D), Matrix{Float64}(D^2), Float64.(vec(y))
end

"""Sample a finite-interval base solution onto a finite or semi-infinite grid.

The thermal/base-flow BVP is integrated through `y=40`.  Grid points beyond
that interval, including the endpoint at infinity, are assigned the governing
far-field state rather than extrapolated from a spline.
"""
function sample_base_state(source, T0, H0, y::Vector{Float64}, p::PhysicalParameters)
    splines = map(values -> BSplineKit.interpolate(
        source.y, vec(values), BSplineOrder(4),
    ), (source.u, source.v, H0, T0))
    source_end = last(source.y)
    F = similar(y)
    G = similar(y)
    H = similar(y)
    T = similar(y)
    for index in eachindex(y)
        point = y[index]
        if isfinite(point) && point <= source_end
            F[index] = splines[1](point)
            G[index] = splines[2](point)
            H[index] = splines[3](point)
            T[index] = splines[4](point)
        else
            F[index] = 0.0
            G[index] = 1.0
            H[index] = last(H0)
            T[index] = 1.0
        end
    end
    # Enforce the prescribed wall and asymptotic states exactly at endpoints.
    F[1], G[1], H[1], T[1] = 0.0, 0.0, 0.0, p.Tw
    F[end], G[end], H[end], T[end] = 0.0, 1.0, last(H0), 1.0
    return F, G, H, T, 1.0 ./ T
end

"""Discretisation and targeted polynomial-eigensolver controls."""
Base.@kwdef struct NumericalParameters
    N::Int = 69
    domain::Symbol = :rational
    ymax::Float64 = 40.0
    target::ComplexF64 = 0.38482 + 0im
    neigs::Int = 2
    tol::Float64 = 1.0e-11
    maxit::Int = 800
    balance::Bool = true
    balance_iterations::Int = 2
    regularization::Float64 = 0.0
    property_perturbations::Bool = true
    base_property_variation::Bool = true
    # :halfline refuses a finite-domain thermal profile when the asymptotic
    # equation cannot satisfy the prescribed wall and far-field temperatures.
    # :finite_diagnostic retains the old calculation with an explicit label.
    thermal_model::Symbol = :halfline
    # Empirical minimum from the registered low-Ro domain sweep: about ten
    # far-field e-folds gives a <1e-3 f-tail fraction for the sampled cases.
    min_thermal_decay_lengths::Float64 = 10.0
end

"""Build a BEK case at fixed system-rotation Mach number.

With Mr=r*abs(DeltaOmega)/a_inf and Ro=DeltaOmega/Omega, the physically
matched path has Mr=abs(Ro)*Msystem. Its exact Ro=0 endpoint has Mr=0 and
requires an incompressible stability solver rather than the singular
compressible pressure scaling `1/Ma^2`.
"""
function at_fixed_system_mach(p::PhysicalParameters, Ro::Real, Msystem::Real)
    isfinite(Ro) && isfinite(Msystem) && Msystem >= 0 ||
        throw(ArgumentError("Ro must be finite and Msystem finite/nonnegative"))
    return PhysicalParameters(
        Tw=p.Tw, Mr=abs(Float64(Ro))*Float64(Msystem), Ro=Float64(Ro),
        R=p.R, beta=p.beta, omega=p.omega, gamma=p.gamma, Pr=p.Pr,
    )
end

struct PreparedCase
    physical::PhysicalParameters
    numerical::NumericalParameters
    Co::Float64
    Ma::Float64
    thermal_bc_error::Float64
    thermal_farfield_exponent::Float64
    thermal_decay_lengths::Float64
    y::Vector{Float64}
    D::Matrix{Float64}
    D2::Matrix{Float64}
    F::Vector{Float64}
    G::Vector{Float64}
    H::Vector{Float64}
    T::Vector{Float64}
    rho::Vector{Float64}
    lambda::Vector{Float64}
    kappa::Vector{Float64}
    evidence_scope::Symbol
end

struct StabilityResult
    case::PreparedCase
    alpha::Vector{ComplexF64}
    vectors::Matrix{ComplexF64}
    backward_errors::Vector{Float64}
    overlaps::Vector{Float64}
end

function check_parameters(p::PhysicalParameters, n::NumericalParameters;
                          allow_unverified_ro::Bool=false)
    p.Tw > 0 || throw(ArgumentError("Tw must be positive"))
    p.Mr >= 0 || throw(ArgumentError("Mr must be nonnegative"))
    p.R > 0 || throw(ArgumentError("R must be positive"))
    p.gamma > 1 || throw(ArgumentError("gamma must exceed one"))
    isapprox(p.Pr, 0.72; atol=0, rtol=8eps(Float64)) || throw(ArgumentError(
        "the current base-flow backend hard-codes Pr=0.72; another Pr requires a rederived base solver",
    ))
    n.N >= 20 || throw(ArgumentError("N must be at least 20"))
    n.domain in (:rational, :finite) || throw(ArgumentError("domain must be :rational or :finite"))
    n.ymax > 0 || throw(ArgumentError("ymax must be positive"))
    n.neigs >= 1 || throw(ArgumentError("neigs must be positive"))
    n.regularization >= 0 || throw(ArgumentError("regularization must be nonnegative"))
    n.thermal_model in (:halfline, :finite_diagnostic) || throw(ArgumentError(
        "thermal_model must be :halfline or :finite_diagnostic"))
    n.min_thermal_decay_lengths >= 0 || throw(ArgumentError(
        "min_thermal_decay_lengths must be nonnegative"))
    if p.Ro != -1.0 && !allow_unverified_ro
        throw(ArgumentError(
            "Ro != -1 invokes a BEK extension not validated by the 2006 paper; " *
            "pass allow_unverified_ro=true only for explicitly labelled diagnostics",
        ))
    end
    return nothing
end

"""Compute velocity and thermal source profiles on `0 <= y <= ymax`.

This replaces the legacy `f_q` callback, whose endpoint expressions used
array indexing instead of evaluating the BVP solution and therefore did not
enforce the far-field conditions.  Both scalar problems are linear Dirichlet
BVPs and are solved directly with a centered tridiagonal discretisation.
"""
function source_profiles(Ro::Float64, Pr::Float64; ymax::Float64=40.0)
    ymax > 0 || throw(ArgumentError("thermal source ymax must be positive"))
    cache_key = (Ro, Pr, ymax)
    return get!(_SOURCE_PROFILE_CACHE, cache_key) do
        base_ymax = 40.0
        source_N = max(10001, round(Int, 250*ymax)+1)
        source_y = range(0.0, ymax; length=source_N)
        returned_u, returned_v, returned_w, returned_du, returned_dv, _ =
            CRD.CRD_BF.sol_baseflowODE(Ro)
        # The preserved backend returns (-U,-V,-W). Undo that historical
        # output reversal once, so every Ro uses Lingwood's direct (U,V,W)
        # convention and no endpoint-specific profile flip remains.
        raw_u, raw_v, raw_w = -returned_u, -returned_v, -returned_w
        du0, dv0 = -returned_du, -returned_dv
        base_y = range(0.0, base_ymax; length=length(raw_u))
        phi = CRD.CRD_BF.phi_var(raw_u, step(base_y), length(base_y))
        base_u, _, base_v, _, base_w, Fu, Fdu, Fdv, Fw, _ = CRD.CRD_BF.velocity(
            raw_u, raw_v, raw_w, du0, dv0, base_y, phi,
        )
        Fv = BSplineKit.interpolate(base_y, base_v, BSplineOrder(4))
        u = [point <= base_ymax ? Fu(point) : 0.0 for point in source_y]
        v = [point <= base_ymax ? Fv(point) : 1.0 for point in source_y]
        w = [point <= base_ymax ? Fw(point) : last(raw_w) for point in source_y]

        function solve_dirichlet(a,b,s,left_value,right_value)
            count = length(source_y)
            interior = count-2
            spacing = step(source_y)
            inverse_h2 = inv(spacing^2)
            lower = zeros(Float64,interior-1)
            diagonal = zeros(Float64,interior)
            upper = zeros(Float64,interior-1)
            rhs = zeros(Float64,interior)
            for row in 1:interior
                index = row+1
                left_coefficient = inverse_h2 + a[index]/(2spacing)
                diagonal[row] = -2inverse_h2-b[index]
                right_coefficient = inverse_h2-a[index]/(2spacing)
                rhs[row] = s[index]
                row > 1 && (lower[row-1] = left_coefficient)
                row < interior && (upper[row] = right_coefficient)
                row == 1 && (rhs[row] -= left_coefficient*left_value)
                row == interior && (rhs[row] -= right_coefficient*right_value)
            end
            values = zeros(Float64,count)
            values[1],values[end] = left_value,right_value
            values[2:end-1] .= Tridiagonal(lower,diagonal,upper) \ rhs
            return values
        end

        sampled_u = u
        sampled_w = w
        sampled_du = [point <= base_ymax ? Fdu(point) : 0.0 for point in source_y]
        sampled_dv = [point <= base_ymax ? Fdv(point) : 0.0 for point in source_y]
        # Unified BEK energy equation in the signed DeltaOmega convention:
        # f'' - Pr*Ro*W*f' - 2Pr*Ro*U*f = 2Pr*(U'^2+V'^2)
        # q'' - Pr*Ro*W*q' = 0.
        # solve_dirichlet represents y'' - transport*y' - reaction*y = rhs.
        transport = Pr*Ro .* sampled_w
        reaction = 2Pr*Ro .* sampled_u
        dissipation = 2Pr .* (sampled_du.^2 .+ sampled_dv.^2)
        f = solve_dirichlet(transport,reaction,dissipation,0.0,0.0)
        q = solve_dirichlet(transport,zeros(source_N),zeros(source_N),1.0,0.0)
        endpoint_error = maximum(abs.((f[1], f[end], q[1] - 1, q[end])))
        endpoint_error <= 1e-8 || error(
            "thermal two-point BVP failed its endpoint Gate: error=$endpoint_error",
        )
        (
            u=Float64.(vec(u)), v=Float64.(vec(v)), w=Float64.(vec(w)),
            f=f, q=q, y=Float64.(collect(source_y)),
            endpoint_error=endpoint_error,
        )
    end
end

"""Build the base state on the paper's similarity coordinate `y`.

This intentionally calls `interp(..., "sim")`; using `"phy"` inserts an
additional coordinate transformation and is inconsistent with equations
(10)--(15) of the paper.
"""
function prepare_case(p::PhysicalParameters=PhysicalParameters(),
                      n::NumericalParameters=NumericalParameters();
                      allow_unverified_ro::Bool=false)
    check_parameters(p, n; allow_unverified_ro=allow_unverified_ro)
    thermal_active = p.Mr > 1e-12 || abs(p.Tw - 1.0) > 1e-12
    thermal_degenerate = abs(p.Ro) <= 1e-12 && thermal_active
    if thermal_degenerate && n.thermal_model == :halfline
        throw(ArgumentError(
            "exact Ro=0 with Mr>0 or Tw!=1 has no prescribed-temperature " *
            "half-line solution: q''=0 and f''=2Pr*(U'^2+V'^2). " *
            "Use thermal_model=:finite_diagnostic only for a labelled " *
            "finite-domain test, or take a matched low-Ro/Mr limit.",
        ))
    end
    Co = 2 - p.Ro - p.Ro^2
    Ma = p.Mr / p.R
    source = source_profiles(p.Ro, p.Pr; ymax=n.ymax)
    # The far-field homogeneous heat equation has exponents 0 and
    # Pr*Ro*W_inf. A nonnegative second exponent cannot provide the
    # decaying thermal mode required by the imposed far-field Dirichlet data.
    farfield_exponent = p.Pr * p.Ro * last(source.w)
    thermal_outflow = !thermal_degenerate && thermal_active &&
        farfield_exponent >= -1e-10
    if thermal_outflow && n.thermal_model == :halfline
        throw(ArgumentError(
            "thermal far field has no decaying similarity mode " *
            "(Ro*W_inf=$(p.Ro*last(source.w))). " *
            "Use thermal_model=:finite_diagnostic only for a labelled " *
            "finite-domain test; a non-similar/outflow thermal model is needed.",
        ))
    end
    decay_lengths = -farfield_exponent*n.ymax
    thermal_domain_unresolved = !thermal_degenerate && !thermal_outflow &&
        thermal_active && decay_lengths < n.min_thermal_decay_lengths
    if thermal_domain_unresolved && n.thermal_model == :halfline
        throw(ArgumentError(
            "thermal domain spans only $decay_lengths far-field decay lengths; " *
            "at least $(n.min_thermal_decay_lengths) are required by this " *
            "numerical Gate. Increase ymax or use " *
            "thermal_model=:finite_diagnostic for a labelled test.",
        ))
    end
    D, D2, y = stability_grid(n)
    T0 = 1 .- ((p.gamma - 1) / 2) .* p.Mr^2 .* source.f .+
         (p.Tw - 1) .* source.q
    H0 = source.w .* T0
    F, G, H, T, rho = sample_base_state(source, T0, H0, y, p)

    Tv = Float64.(T)
    # At exact Ro=0 the thermal similarity equations lose their normal
    # transport. With nontrivial Mr or Tw this finite-domain construction is
    # a singular diagnostic, not a half-line base state.
    evidence_scope = thermal_degenerate ? :ro0_thermal_degenerate :
        thermal_outflow ? :thermal_outflow_unverified :
        thermal_domain_unresolved ? :thermal_domain_unresolved :
        :derived_bek_candidate
    return PreparedCase(
        p, n, Co, Ma, source.endpoint_error,
        farfield_exponent, decay_lengths,
        y, D, D2, F, G, H, Tv, rho,
        -(2 / 3) .* Tv, Tv ./ p.Pr, evidence_scope,
    )
end

"""Build coefficient tensors in Lingwood's direct `(U,V,W)` convention."""
function corrected_coefficients(c::PreparedCase)
    p, n = c.physical, c.numerical
    coefficient_at = ro -> CRD.Spatial_mode_BEK(
        c.F, c.G, c.H, c.rho, c.lambda, c.kappa, c.T,
        p.Pr, p.gamma, p.R, c.Ma, n.N, ro, c.Co, c.D, c.D2;
        regularization=n.regularization,
        property_perturbations=n.property_perturbations,
        base_property_variation=n.base_property_variation,
    )
    coefficient_zero = coefficient_at(0.0)
    coefficient_one = coefficient_at(1.0)
    components = map(fieldnames(CRD.COF)) do name
        z = getfield(coefficient_zero, name)
        z + p.Ro * (getfield(coefficient_one, name) - z)
    end
    coefficients = CRD.COF(components...)

    # General-BEK corrections that cannot be obtained by inserting Ro and Co
    # into the published rotating-disk operator. Field ordering is
    # (rho,u,v,w,T). In direct Lingwood variables their coefficient is Ro.
    npoints = n.N + 1
    rr = 1:npoints
    ru = npoints+1:2npoints
    rT = 4npoints+1:5npoints
    eye = Matrix{Float64}(I, npoints, npoints)
    D1 = copy(coefficients.D1)
    C = copy(coefficients.C)
    dC = copy(coefficients.dC)

    # Compressible continuity: cylindrical u/r and base W transport carry Ro
    # in the positive-R local BEK scaling.
    continuation_factor = p.Ro
    D1[rr,ru] .= continuation_factor .* c.rho .* eye
    D1[rr,rr] .= continuation_factor .* (2 .* c.F .+ c.rho .* (c.D*c.H)) .* eye
    C[rr,rr] .= continuation_factor .* c.rho .* c.H .* eye
    dC[rr,rr] .= c.D * Diagonal(diag(C[rr,rr]))

    # Energy: the base wall-normal material derivative and its pressure-work
    # partners inherit the same Ro as Ro*W*D in the momentum equations.
    property_scale = n.property_perturbations ? 1.0 : 0.0
    Drho = c.D*c.rho
    DT = c.D*c.T
    D2T = c.D2*c.T
    DU = c.D*c.F
    DV = c.D*c.G
    C[rT,rr] .= -(p.gamma-1)/p.gamma .* continuation_factor .* c.rho .* c.H .* c.T .* eye
    C[rT,rT] .= (
        (continuation_factor/p.gamma) .* c.rho.^2 .* c.H .+
        (property_scale/p.Pr) .* c.rho.^2 .* DT
    ) .* eye
    dC[rT,rr] .= c.D * Diagonal(diag(C[rT,rr]))
    dC[rT,rT] .= c.D * Diagonal(diag(C[rT,rT]))
    D1[rT,rr] .= (continuation_factor/p.gamma) .* c.rho .* c.H .* DT .* eye
    D1[rT,rT] .= (
        -(p.gamma-1)/p.gamma .* continuation_factor .* c.rho .* c.H .* Drho .-
        (property_scale/p.Pr) .* (c.rho .* Drho .* DT .+ c.rho.^2 .* D2T) .-
        property_scale .* (p.gamma-1) .* c.Ma^2 .* c.rho.^2 .* (DU.^2 .+ DV.^2)
    ) .* eye

    corrected_components = map(fieldnames(CRD.COF)) do name
        name === :D1 ? D1 : name === :C ? C : name === :dC ? dC :
            getfield(coefficients, name)
    end
    coefficients = CRD.COF(corrected_components...)
    return coefficients
end

"""Assemble the boundary-reduced quadratic eigenproblem
`(A0 + alpha*A1 + alpha^2*A2) q = 0`."""
function assemble_qep(c::PreparedCase)
    c.Ma > 0 || throw(ArgumentError(
        "Mr=0 makes the compressible pressure scaling 1/Ma^2 singular; " *
        "use the incompressible BEK operator at exact Ro=0, or approach " *
        "it with Mr=abs(Ro)*Msystem at nonzero Ro.",
    ))
    p, n = c.physical, c.numerical
    coefficients = corrected_coefficients(c)
    raw = CRD.assemble_mat(coefficients, c.D, c.D2, p.beta, p.omega)
    return CRD.boundary_condition(raw..., n.N)
end

"""Return exact residuals for the explicit Ro/Co placement gates."""
function operator_structure_gates(c::PreparedCase)
    p, n = c.physical, c.numerical
    cof = corrected_coefficients(c)
    np = n.N+1
    rr, ru, rv, rT = 1:np, np+1:2np, 2np+1:3np, 4np+1:5np
    gamma_block = 2p.Ro .* c.G .+ c.Co
    eye = Matrix{Float64}(I,np,np)
    property_scale = n.property_perturbations ? 1.0 : 0.0
    continuation_factor = p.Ro
    expected_energy = (
        (continuation_factor/p.gamma) .* c.rho.^2 .* c.H .+
        (property_scale/p.Pr) .* c.rho.^2 .* (c.D*c.T)
    ) .* eye
    return (
        continuity_curvature=norm(cof.D1[rr,ru] .-
                                  continuation_factor .* c.rho .* eye),
        continuity_base_transport=norm(
            cof.C[rr,rr] .- continuation_factor .* c.rho .* c.H .* eye),
        radial_coriolis=norm(diag(cof.D1[ru,rv]) .+
                             c.rho .* gamma_block),
        azimuthal_coriolis=norm(diag(cof.D1[rv,ru]) .-
                                c.rho .* gamma_block),
        bodewadt_direct_Co=abs(p.Ro-1) < 10eps() ? abs(c.Co) : NaN,
        ekman_curvature_should_vanish=abs(p.Ro) < 10eps() ?
            norm(cof.D1[rr,ru]) : NaN,
        energy_wallnormal_scale=norm(cof.C[rT,rT] .- expected_energy),
    )
end

"""Solve for modes nearest `numerical.target` using a balanced targeted QEP."""
function solve_mode(c::PreparedCase; initial_vector=nothing)
    n = c.numerical
    A0, A1, A2 = assemble_qep(c)
    if n.balance
        B0, B1, B2, _, right_scale = CRD.balance_polynomial(
            A0, A1, A2; iterations=n.balance_iterations,
        )
    else
        B0, B1, B2 = ComplexF64.(A0), ComplexF64.(A1), ComplexF64.(A2)
        right_scale = ones(Float64, size(A0, 1))
    end
    start_vector = isnothing(initial_vector) ? ones(ComplexF64, size(A0, 1)) :
        ComplexF64.(vec(initial_vector)) ./ right_scale
    values, balanced_vectors = CRD.iar(
        CRD.PEP([B0, B1, B2]); σ=n.target, neigs=n.neigs,
        maxit=n.maxit, tol=n.tol, v=start_vector,
    )
    vectors = reshape(right_scale, :, 1) .* balanced_vectors
    for column in axes(vectors, 2)
        vectors[:, column] ./= norm(view(vectors, :, column))
    end
    overlaps = isnothing(initial_vector) ? fill(NaN, length(values)) :
        [CRD.mode_overlap(initial_vector, view(vectors, :, j)) for j in eachindex(values)]
    order = isnothing(initial_vector) ?
        sortperm(eachindex(values); by=j -> abs(values[j] - n.target)) :
        sortperm(eachindex(values); by=j -> (-overlaps[j], abs(values[j] - n.target)))
    values, vectors, overlaps = values[order], vectors[:, order], overlaps[order]
    residuals = [
        CRD.polynomial_backward_error(A0, A1, A2, values[j], vectors[:, j])
        for j in eachindex(values)
    ]
    return StabilityResult(
        c, ComplexF64.(values), ComplexF64.(vectors),
        Float64.(residuals), Float64.(overlaps),
    )
end

"""Numerical and physical consistency diagnostics for a solved case."""
function validate_result(result::StabilityResult)
    c = result.case
    A0, A1, A2 = assemble_qep(c)
    recomputed = [
        CRD.polynomial_backward_error(A0, A1, A2, result.alpha[j], result.vectors[:, j])
        for j in eachindex(result.alpha)
    ]
    base_bc = (
        wall_F=abs(c.F[1]), wall_G=abs(c.G[1]), wall_H=abs(c.H[1]),
        wall_T=abs(c.T[1] - c.physical.Tw),
        far_F=abs(c.F[end]), far_G=abs(c.G[end] - 1), far_T=abs(c.T[end] - 1),
    )
    finite_modes = all(isfinite, real.(result.alpha)) && all(isfinite, imag.(result.alpha))
    return (
        finite_modes=finite_modes,
        base_bc=base_bc,
        stored_backward_errors=result.backward_errors,
        recomputed_backward_errors=recomputed,
        matrix_size=size(A0),
        Ma_identity_error=abs(c.Ma * c.physical.R - c.physical.Mr),
        thermal_bc_error=c.thermal_bc_error,
        evidence_scope=c.evidence_scope,
    )
end

"""Retired legacy compact/expanded comparison.

The expanded routine contains the endpoint-specific `Ro=-1 -> +1` rewrite
that the direct-variable BEK operator removes.  Comparing it to the corrected
operator would therefore test two different variable conventions rather than
the same equations.
"""
function audit_expanded_operator(::PreparedCase)
    throw(ArgumentError(
        "audit_expanded_operator is retired for the direct-variable BEK operator; " *
        "use operator_structure_gates and directional_qep_test",
    ))
end

"""Independent matrix-polynomial directional-derivative check.

For a fixed test vector `q`, compares the analytic derivative
`(A1 + 2*alpha*A2)q` with a centered finite difference of
`P(alpha)q`.  This checks the assembled QEP algebra separately from the
eigensolver residual and is useful when changing the general-Ro blocks.
"""
function directional_qep_test(c::PreparedCase;
                              alpha::ComplexF64=ComplexF64(c.numerical.target),
                              h::Float64=1e-6,
                              seed::Int=20260911)
    h > 0 || throw(ArgumentError("h must be positive"))
    A0, A1, A2 = assemble_qep(c)
    rng = Random.MersenneTwister(seed)
    q = randn(rng, ComplexF64, size(A0, 1))
    P(a) = (A0 + a .* A1 + (a^2) .* A2) * q
    analytic = (A1 + 2alpha .* A2) * q
    finite_difference = (P(alpha + h) - P(alpha - h)) / (2h)
    relative_error = norm(finite_difference - analytic) /
        max(norm(analytic), eps(Float64))
    return (relative_error=relative_error, h=h, alpha=alpha,
            vector_norm=norm(q))
end

function write_result(result::StabilityResult, output_dir::AbstractString)
    mkpath(output_dir)
    c, p, n = result.case, result.case.physical, result.case.numerical
    writedlm(joinpath(output_dir, "base_profile.csv"),
             hcat(c.y, c.F, c.G, c.H, c.T, c.rho), ',')
    open(joinpath(output_dir, "eigenvalues.csv"), "w") do io
        println(io, "mode,alpha_real,alpha_imag,backward_error,overlap")
        for j in eachindex(result.alpha)
            println(io, join((j, real(result.alpha[j]), imag(result.alpha[j]),
                              result.backward_errors[j], result.overlaps[j]), ','))
        end
    end
    open(joinpath(output_dir, "metadata.txt"), "w") do io
        println(io, "model=compressible_BEK_spatial_LSA")
        println(io, "evidence_scope=$(c.evidence_scope)")
        println(io, "thermal_farfield_exponent=$(c.thermal_farfield_exponent)")
        println(io, "thermal_decay_lengths=$(c.thermal_decay_lengths)")
        println(io, "coordinate=Dorodnitsyn_Howarth_similarity_y")
        for name in fieldnames(PhysicalParameters)
            println(io, "$(name)=$(getfield(p, name))")
        end
        println(io, "Ma_internal=$(c.Ma)")
        println(io, "Mr_definition=local_differential_rotation_Mach_r_absDeltaOmega_over_ainf")
        println(io, "Co=$(c.Co)")
        println(io, "thermal_bc_error=$(c.thermal_bc_error)")
        for name in fieldnames(NumericalParameters)
            println(io, "$(name)=$(getfield(n, name))")
        end
    end
    return output_dir
end

end

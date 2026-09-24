module TraditionalBoussinesqTemporal

using LinearAlgebra

export temporal_operator, leading_temporal_eigenvalue

"""Chebyshev antiderivative with every primitive equal to zero at z=0."""
function integration_matrix(z::AbstractVector)
    n = length(z)
    x = 2 .* z .- 1
    theta = acos.(clamp.(x, -1.0, 1.0))
    V = hcat((cos.(degree .* theta) for degree in 0:n-1)...)
    Iint = zeros(n, n)
    Iint[:, 1] .= 0.5 .* (x .+ 1)
    n >= 2 && (Iint[:, 2] .= 0.25 .* (x.^2 .- 1))
    for degree in 2:n-1
        tplus = cos.((degree + 1) .* theta)
        tminus = cos.((degree - 1) .* theta)
        at_left = 0.5 * (((-1.0)^(degree + 1)) / (degree + 1) -
                          ((-1.0)^(degree - 1)) / (degree - 1))
        Iint[:, degree + 1] .= 0.5 .* (
            0.5 .* (tplus ./ (degree + 1) .- tminus ./ (degree - 1)) .-
            at_left)
    end
    return Iint / V
end

function linear_interpolate(x, y, targets)
    out = similar(targets, Float64)
    for (i, target) in pairs(targets)
        k = clamp(searchsortedlast(x, target), 1, length(x) - 1)
        a = (target - x[k]) / (x[k + 1] - x[k])
        out[i] = (1 - a) * y[k] + a * y[k + 1]
    end
    return out
end

"""
Temporal operator for axisymmetric, similarity-preserving disturbances.

The perturbation state is (f, g, theta); continuity determines h and the
zero-net-radial-flow constraint removes the pressure direction. This is the
temporal stability problem used for the saddle-node branch classification,
not the full three-dimensional cavity spectrum.
"""
function temporal_operator(profile, degree::Int, re_h::Real, pr::Real;
                           closure::Symbol=:traditional)
    closure in (:traditional, :blackburn) ||
        throw(ArgumentError("closure must be :traditional or :blackburn"))
    n = degree + 1
    z = 0.5 .* (1 .- cos.(pi .* collect(0:degree) ./ degree))
    c = vcat(2.0, ones(degree - 1), 2.0) .* ((-1.0) .^ collect(0:degree))
    X = repeat(cos.(pi .* collect(0:degree) ./ degree), 1, n)
    dX = X .- X'
    D = (c * (1.0 ./ c)') ./ (dX + Matrix{Float64}(I, n, n))
    D .-= Diagonal(vec(sum(D; dims=2)))
    D = -2.0 .* D
    D2 = D * D
    integration = integration_matrix(z)
    weights = vec(integration[end, :])
    values = profile.fields
    fields = [linear_interpolate(profile.z, vec(values[row, :]), z) for row in 1:4]
    h, f, g, temperature = fields
    fz, gz, tz = D * f, D * g, D * temperature
    interior = 2:n-1
    m = length(interior)
    extension = zeros(n, m)
    extension[interior, :] .= Matrix{Float64}(I, m, m)
    weights_i = weights[interior]
    q = nullspace(reshape(weights_i, 1, :))
    size(q) == (m, m - 1) || error("unexpected zero-integral basis dimension")
    pressure_projection = Matrix{Float64}(I, m, m) .-
                          ones(m) * transpose(weights_i ./ sum(weights_i))
    f_map = extension * q
    g_map = extension
    temperature_map = extension
    h_map = -2sqrt(re_h) .* integration * f_map
    f_i, g_i, h_i = f[interior], g[interior], h[interior]
    fz_i, gz_i, tz_i = fz[interior], gz[interior], tz[interior]
    d1_f, d2_f = D * f_map, D2 * f_map
    d1_g, d2_g = D * g_map, D2 * g_map
    d1_t, d2_t = D * temperature_map, D2 * temperature_map
    f_i_map, h_i_map = f_map[interior, :], h_map[interior, :]
    sqrt_re = sqrt(re_h)
    if closure == :traditional
        radial_f = d2_f[interior, :] ./ re_h .-
                   h_i .* d1_f[interior, :] ./ sqrt_re .-
                   fz_i .* h_i_map ./ sqrt_re .- 2f_i .* f_i_map
        radial_g = 2 .* Diagonal(g_i)
        radial_t = -Matrix{Float64}(I, m, m)
        azimuthal_f = -gz_i .* h_i_map ./ sqrt_re .- 2g_i .* f_i_map
        azimuthal_g = d2_g[interior, :] ./ re_h .-
                      h_i .* d1_g[interior, :] ./ sqrt_re .- 2 .* Diagonal(f_i)
        azimuthal_t = zeros(m, m)
        mass_f = Matrix{Float64}(I, m - 1, m - 1)
        mass_g = Matrix{Float64}(I, m, m)
    else
        chi_i = 2.0 .- temperature[interior]
        # Blackburn weights the complete material acceleration by chi=2-T.
        # Perturbing chi supplies the temperature blocks through the base
        # radial and azimuthal accelerations.
        ar = h_i .* fz_i ./ sqrt_re .+ f_i.^2 .- g_i.^2 .+ profile.pressure
        atheta = h_i .* gz_i ./ sqrt_re .+ 2 .* f_i .* g_i
        radial_f = d2_f[interior, :] ./ re_h .-
                   chi_i .* h_i .* d1_f[interior, :] ./ sqrt_re .-
                   chi_i .* fz_i .* h_i_map ./ sqrt_re .-
                   2 .* chi_i .* f_i .* f_i_map
        radial_g = 2 .* Diagonal(chi_i .* g_i)
        radial_t = Diagonal(ar)
        azimuthal_f = -chi_i .* gz_i .* h_i_map ./ sqrt_re .-
                      2 .* chi_i .* g_i .* f_i_map
        azimuthal_g = d2_g[interior, :] ./ re_h .-
                      chi_i .* h_i .* d1_g[interior, :] ./ sqrt_re .-
                      2 .* chi_i .* f_i .* g_map[interior, :]
        azimuthal_t = Diagonal(atheta)
        mass_f = transpose(q) * Diagonal(chi_i) * q
        mass_g = Diagonal(chi_i)
    end
    thermal_f = -tz_i .* h_i_map ./ sqrt_re
    thermal_t = d2_t[interior, :] ./ (re_h * pr) .-
                h_i .* d1_t[interior, :] ./ sqrt_re
    radial_f_reduced = transpose(q) * pressure_projection * radial_f
    radial_g_reduced = transpose(q) * pressure_projection * radial_g
    radial_t_reduced = transpose(q) * pressure_projection * radial_t
    operator = [radial_f_reduced radial_g_reduced radial_t_reduced;
                azimuthal_f azimuthal_g azimuthal_t;
                thermal_f zeros(m, m) thermal_t]
    mass = [mass_f zeros(m - 1, m) zeros(m - 1, m);
            zeros(m, m - 1) mass_g zeros(m, m);
            zeros(m, m - 1) zeros(m, m) Matrix{Float64}(I, m, m)]
    spectrum = closure == :traditional ? eigen(operator) : eigen(operator, mass)
    order = sortperm(eachindex(spectrum.values);
                     by=i -> (-real(spectrum.values[i]), abs(imag(spectrum.values[i]))) )
    return (; operator, eigenvalues=ComplexF64.(spectrum.values[order]),
            eigenvectors=ComplexF64.(spectrum.vectors[:, order]), z,
            f_map, h_map, g_map, temperature_map)
end

leading_temporal_eigenvalue(result) = first(result.eigenvalues)

end

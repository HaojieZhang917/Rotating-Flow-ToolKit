#!/usr/bin/env julia

using Printf
using Test
using LinearAlgebra
using Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))
using CompressRotatingCavity
include(joinpath(ROOT, "src", "CRC_RA.jl"))
include(joinpath(ROOT, "src", "CRC_MEAN.jl"))
using .CRC_RA
using .CRC_MEAN

const N = 65
const RE = 2000.0
const R = 330.0
const AS = 0.0
const nmode = 27
const beta = nmode / R
const omega = 0.0
const kappa_max = 4.0
const Nk = 41
const ls = 0.5
const h0 = 1.0 / RE
const Mtheta = 4nmode
const theta = 2pi .* collect(0:Mtheta-1) ./ Mtheta
const kappa_min = 0.05
const ds = (kappa_max - kappa_min) / (Nk - 1)
const kappas = collect(range(kappa_min, kappa_max, length=Nk))

poly(L0, L1, L2, a) = L0 + a * L1 + a^2 * L2
gaussian_hat(k; h) = h * sqrt(pi / ls) * exp(-k^2 / (4ls))

function boundary_forcing(flow, N, k; height)
    hhat = gaussian_hat(k; h=height)
    return wall_boundary_vector(N;
        uw=-(flow.D * flow.F)[1] * hhat,
        vw=-(flow.D * flow.G)[1] * hhat)
end

function trap_weight(j, n)
    return (j == 1 || j == n) ? 0.5 : 1.0
end

function reconstruct_mode_fields(qhats, dzhats, hhat, D, s, h)
    n = size(dzhats, 1) ÷ 3
    ur = zeros(Float64, n, Mtheta); ut = similar(ur); uz = similar(ur)
    dr = zeros(Float64, n, Mtheta); dt = similar(dr); dz = similar(dr)
    hp = zeros(Float64, Mtheta)
    for j in eachindex(kappas)
        k = kappas[j]
        weight = ds * trap_weight(j, length(kappas)) / pi
        phase = exp.(im .* (k*s .+ nmode .* theta))
        q = qhats[:, j]
        dq = dzhats[:, j]
        ur .+= weight .* real.(q[1:n] .* transpose(phase))
        ut .+= weight .* real.(q[n+1:2n] .* transpose(phase))
        uz .+= weight .* real.(q[2n+1:3n] .* transpose(phase))
        dr .+= weight .* real.((im*k .* q[1:n]) .* transpose(phase))
        dt .+= weight .* real.((im*nmode/R .* q[1:n]) .* transpose(phase))
        dz .+= weight .* real.(dq[1:n] .* transpose(phase))
        # Wall displacement reconstructed with the same Fourier convention.
        hp .+= weight .* real.(hhat[j] .* phase)
    end
    return (; ur, ut, uz, dr, dt, dz, hp)
end

function mean_acceleration(fields, D, R)
    ur, ut, uz = fields.ur, fields.ut, fields.uz
    dr = fields.dr
    dt = fields.dt
    dzur = D * ur
    dzut = D * ut
    dzuz = D * uz
    ar = mean(ur .* dr .+ ut ./ R .* dt .+ uz .* dzur .- ut.^2 ./ R; dims=2)[:, 1]
    atheta = mean(ur .* (fields.dr .* 0 .+ 0) .+ 0; dims=2)[:, 1] # overwritten below
    # Reconstruct azimuthal/radial derivatives of the separate components.
    # The dt field above is d(theta)ur; the corresponding component fields are
    # generated from the Fourier data in the caller and stored separately.
    return ar, atheta, mean(ur .* 0; dims=2)[:, 1]
end

println("Full radial-spectrum stationary roughness forcing verification")
println("  Re_h=$RE, R=$R, n=$nmode, beta=$beta, omega=$omega, Nk=$Nk, kmax=$kappa_max")

flow = (; baseflow(RE, -1.0, AS, 1, N)..., Res=RE)
solution0 = normalize_modes(eigsolve(flow, RE, -1.0, 1, 0.4, R, omega, beta))
alpha0 = solution0.eigval[1]
sta = CompressRotatingCavity.CRC_STA
cof = sta.Spatial_mode_BEK1(flow.F, flow.G .- 1, flow.H, R, N,
                             flow.D, flow.D2, RE)
L0, L1, L2 = sta.assemble_mat(cof, flow.D, flow.D2, beta, omega, R)

qhats = zeros(ComplexF64, 3(N+1), length(kappas))
dzhats = zeros(ComplexF64, 3(N+1), length(kappas))
hhat = [gaussian_hat(k; h=h0) for k in kappas]
for (j, k) in enumerate(kappas)
    P = poly(L0, L1, L2, k)
    q = lifted_response(P, N, boundary_forcing(flow, N, k; height=h0))
    all(isfinite, q) || println("nonfinite qhat at kappa=", k)
    qhats[:, j] .= q[1:3(N+1)]
    dzhats[:, j] .= vcat(flow.D*q[1:N+1],
                         flow.D*q[N+2:2(N+1)],
                         flow.D*q[2(N+1)+1:3(N+1)])
    all(isfinite, dzhats[:, j]) || println("nonfinite dzhat at kappa=", k)
end

function fields_at(s)
    n = N + 1
    ur = zeros(Float64, n, Mtheta); ut = zeros(Float64, n, Mtheta); uz = zeros(Float64, n, Mtheta)
    drur = zeros(Float64, n, Mtheta); drut = zeros(Float64, n, Mtheta); druz = zeros(Float64, n, Mtheta)
    dthur = zeros(Float64, n, Mtheta); dthut = zeros(Float64, n, Mtheta); dthuz = zeros(Float64, n, Mtheta)
    dzur = zeros(Float64, n, Mtheta); dzut = zeros(Float64, n, Mtheta); dzuz = zeros(Float64, n, Mtheta)
    hphys = zeros(Float64, Mtheta)
    for j in eachindex(kappas)
        k = kappas[j]
        weight = ds * trap_weight(j, length(kappas)) / pi
        phase = exp.(im .* (k*s .+ nmode .* theta))
        q = qhats[:, j]; dq = dzhats[:, j]
        for (arr, qblock) in ((ur,q[1:n]), (ut,q[n+1:2n]), (uz,q[2n+1:3n]))
            arr .+= weight .* real.(qblock .* transpose(phase))
        end
        drur .+= weight .* real.((im*k .* q[1:n]) .* transpose(phase))
        drut .+= weight .* real.((im*k .* q[n+1:2n]) .* transpose(phase))
        druz .+= weight .* real.((im*k .* q[2n+1:3n]) .* transpose(phase))
        dthur .+= weight .* real.((im*nmode/R .* q[1:n]) .* transpose(phase))
        dthut .+= weight .* real.((im*nmode/R .* q[n+1:2n]) .* transpose(phase))
        dthuz .+= weight .* real.((im*nmode/R .* q[2n+1:3n]) .* transpose(phase))
        dzur .+= weight .* real.(dq[1:n] .* transpose(phase))
        dzut .+= weight .* real.(dq[n+1:2n] .* transpose(phase))
        dzuz .+= weight .* real.(dq[2n+1:3n] .* transpose(phase))
        hphys .+= weight .* real.(hhat[j] .* phase)
    end
    ar = mean(ur .* drur .+ ut ./ R .* dthur .+ uz .* dzur .- ut.^2 ./ R; dims=2)[:,1]
    at = mean(ur .* drut .+ ut ./ R .* dthut .+ uz .* dzut .+ ur .* ut ./ R; dims=2)[:,1]
    az = mean(ur .* druz .+ ut ./ R .* dthuz .+ uz .* dzuz; dims=2)[:,1]
    return (; ar, at, az, hphys, ur, ut, uz, dzur, dzut, dzuz)
end

center = fields_at(0.0)
println(@sprintf("center forcing norms: radial=%.6e theta=%.6e axial=%.6e",
    norm(center.ar)/sqrt(N+1), norm(center.at)/sqrt(N+1), norm(center.az)/sqrt(N+1)))
println(@sprintf("center reconstructed wall h: max=%.6e, target=%.6e",
    maximum(abs, center.hphys), h0))

wall2 = (
    -mean(center.hphys .* center.dzur[1, :]) -
        0.5*mean(center.hphys.^2)*(flow.D2*vec(flow.F))[1],
    -mean(center.hphys .* center.dzut[1, :]) -
        0.5*mean(center.hphys.^2)*(flow.D2*vec(flow.G))[1],
    -mean(center.hphys .* center.dzuz[1, :]) -
        0.5*mean(center.hphys.^2)*(flow.D2*vec(flow.H))[1])
corr = mean_similarity_correction(flow, flow.D, RE, center.ar, center.at;
                                   wall_values=wall2)
println(@sprintf("frozen-R integrated correction norm=%.6e, linear solve residual=%.3e",
    sqrt(sum(abs2,corr.F2)+sum(abs2,corr.G2)+sum(abs2,corr.H2))/sqrt(N+1),
    corr.residual))

flow_corr = (; flow..., F=flow.F .+ corr.F2,
             G=flow.G .+ corr.G2, H=flow.H .+ corr.H2)
solution_corr = normalize_modes(eigsolve(flow_corr, RE, -1.0, 1,
                                         real(alpha0), R, omega, beta))
alpha_corr = solution_corr.eigval[1]
Q0 = q_receptivity(solution0, cheb_quad(N))[1]
Qc = q_receptivity(solution_corr, cheb_quad(N))[1]
Cr0 = cr_receptivity(solution0, flow, RE, R, h0, ls, [Q0]).Cr
Crc = cr_receptivity(solution_corr, flow_corr, RE, R, h0, ls, [Qc]).Cr
println(@sprintf("frozen-R spectral alpha: %.12f %+.12fi -> %.12f %+.12fi",
    real(alpha0), imag(alpha0), real(alpha_corr), imag(alpha_corr)))
println(@sprintf("delta alpha=%.6e %+.6ei, delta Cr=%.6e",
    real(alpha_corr-alpha0), imag(alpha_corr-alpha0), Crc-Cr0))

rows = ["z,mean_ar,mean_at,mean_az,F2,G2,H2"]
for j in eachindex(flow.z)
    push!(rows, @sprintf("%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e",
        flow.z[j], center.ar[j], center.at[j], center.az[j],
        corr.F2[j], corr.G2[j], corr.H2[j]))
end
out = joinpath(ROOT, "benchmark_results", "mean_forcing_spectrum_stationary.csv")
write(out, join(rows, "\n") * "\n")
println("Saved: $out")

@test maximum(abs, center.hphys) ≈ h0 rtol=5e-2
@test corr.residual < 1e-8
println("PASS: radial-spectrum forcing and frozen-R correction computed")

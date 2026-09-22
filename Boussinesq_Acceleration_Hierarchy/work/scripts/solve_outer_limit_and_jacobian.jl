#!/usr/bin/env julia

using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTDIR = joinpath(ROOT, "work", "results", "corrected_energy_outer_limit")
const PR = 0.72
const GAMMA = 1.0

function coefficients(r)
    s = (2 - r - r^2) / 2
    omega = s + r
    return s, omega
end

function theta_curve(r; gamma=GAMMA)
    s, omega = coefficients(r)
    1 + (1 - (omega / s)^2) / gamma
end

function theta_curve_derivative(r; gamma=GAMMA)
    s, omega = coefficients(r)
    sp = -0.5 - r
    omegap = sp + 1
    ratio = omega / s
    ratio_prime = (omegap * s - omega * sp) / s^2
    -2ratio * ratio_prime / gamma
end

function rhs(state, r)
    g, p, h = state
    s, omega = coefficients(r)
    q = s + r * g
    chi = omega^2 / q^2
    pp = 3r * p^2 / q - PR * r * h * p
    hp = (pp + chi * r * h * p) / (chi * q)
    [p, pp, hp]
end

function rk4_endpoint(p0, r, ymax, step; save_profile=false)
    nsteps = ceil(Int, ymax / step)
    dy = ymax / nsteps
    state = [0.0, p0, 0.0]
    if save_profile
        profile = Matrix{Float64}(undef, nsteps + 1, 6)
        s, omega = coefficients(r)
        q = s + r * state[1]
        theta = 1 + (1 - omega^2 / q^2) / GAMMA
        profile[1, :] .= (0.0, state[1], state[2], state[3], theta, omega^2 / q^2)
    end
    for k in 1:nsteps
        k1 = rhs(state, r)
        k2 = rhs(state .+ 0.5dy .* k1, r)
        k3 = rhs(state .+ 0.5dy .* k2, r)
        k4 = rhs(state .+ dy .* k3, r)
        state .+= (dy / 6) .* (k1 .+ 2k2 .+ 2k3 .+ k4)
        all(isfinite, state) || error("outer shooting integration diverged")
        if save_profile
            s, omega = coefficients(r)
            q = s + r * state[1]
            theta = 1 + (1 - omega^2 / q^2) / GAMMA
            profile[k + 1, :] .= (k * dy, state[1], state[2], state[3], theta,
                                  omega^2 / q^2)
        end
    end
    save_profile ? (state, profile) : state
end

function shooting_residual(z, ymax, step)
    state = rk4_endpoint(z[1], z[2], ymax, step)
    [state[1] - 1, state[3] + 1]
end

function shooting_jacobian(z, ymax, step; relative_step=2e-6)
    jac = zeros(2, 2)
    for j in 1:2
        delta = relative_step * max(abs(z[j]), 0.1)
        zp = copy(z); zm = copy(z)
        zp[j] += delta; zm[j] -= delta
        jac[:, j] .= (shooting_residual(zp, ymax, step) .-
                      shooting_residual(zm, ymax, step)) ./ (2delta)
    end
    jac
end

function solve_shooting(initial, ymax, step; tolerance=2e-11, max_iterations=20)
    z = copy(Float64.(initial))
    for iteration in 1:max_iterations
        residual = shooting_residual(z, ymax, step)
        norm0 = norm(residual, Inf)
        norm0 < tolerance && return z, residual, iteration
        jac = shooting_jacobian(z, ymax, step)
        dz = -(jac \ residual)
        damping = 1.0
        accepted = false
        while damping >= 2.0^-16
            trial = z .+ damping .* dz
            trial_residual = shooting_residual(trial, ymax, step)
            if all(isfinite, trial_residual) && norm(trial_residual, Inf) < norm0
                z = trial
                accepted = true
                break
            end
            damping *= 0.5
        end
        accepted || error("outer shooting Newton failed at residual $norm0")
    end
    error("outer shooting did not converge")
end

function main()
mkpath(OUTDIR)
configs = [(ymax=20.0, step=0.01),
           (ymax=30.0, step=0.01),
           (ymax=40.0, step=0.005),
           (ymax=60.0, step=0.005),
           (ymax=60.0, step=0.0025)]
initial = [0.75, -0.48]
rows = NamedTuple[]
for config in configs
    z, residual, iterations = solve_shooting(initial, config.ymax, config.step)
    jac = shooting_jacobian(z, config.ymax, config.step)
    singular = svdvals(jac)
    push!(rows, (; config.ymax, config.step, p0=z[1], r=z[2],
                  tw=theta_curve(z[2]), residual=norm(residual, Inf),
                  sigma_max=singular[1], sigma_min=singular[end],
                  condition=singular[1] / singular[end], determinant=det(jac),
                  dTdr=theta_curve_derivative(z[2]), iterations))
    initial = z
    @printf("Ymax=%4.0f dY=%7.4f p0=%+.12f r=%+.12f Tw=%.12f residual=%.2e cond=%.3e\n",
            config.ymax, config.step, z[1], z[2], theta_curve(z[2]),
            norm(residual, Inf), singular[1] / singular[end])
end

open(joinpath(OUTDIR, "convergence_and_jacobian.csv"), "w") do io
    println(io, "Ymax,step,p0,Ro_t,Tw_t,residual,sigma_max,sigma_min,condition,determinant,dTcurve_dRo,iterations")
    for row in rows
        @printf(io, "%.8f,%.8f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%d\n",
                row.ymax, row.step, row.p0, row.r, row.tw, row.residual,
                row.sigma_max, row.sigma_min, row.condition, row.determinant,
                row.dTdr, row.iterations)
    end
end

best = rows[end]
state, profile = rk4_endpoint(best.p0, best.r, best.ymax, best.step; save_profile=true)
open(joinpath(OUTDIR, "outer_profile.csv"), "w") do io
    println(io, "Y,g,g_Y,h,theta,chi")
    for i in axes(profile, 1)
        @printf(io, "%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n", profile[i, :]...)
    end
end

open(joinpath(OUTDIR, "jacobian_fd_sensitivity.csv"), "w") do io
    println(io, "relative_step,sigma_max,sigma_min,condition,determinant")
    zbest = [best.p0, best.r]
    for relative_step in (1e-4, 1e-5, 2e-6, 5e-7, 1e-7)
        jac = shooting_jacobian(zbest, best.ymax, best.step; relative_step)
        singular = svdvals(jac)
        @printf(io, "%.12e,%.12e,%.12e,%.12e,%.12e\n", relative_step,
                singular[1], singular[end], singular[1] / singular[end], det(jac))
    end
end

open(joinpath(OUTDIR, "scope.txt"), "w") do io
    println(io, "Corrected-energy leading outer BVP only; no finite-epsilon base-flow or fold continuation solve.")
    println(io, "Energy equation: theta_YY + Pr*Ro_t*h*theta_Y = 0.")
    println(io, "The reported 2x2 Jacobian is the reduced augmented shooting Jacobian")
    println(io, "d(g(Ymax)-1,h(Ymax)+1)/d(g_Y(0),Ro), after field IVP elimination.")
    println(io, "It tests isolation of the normalized outer BVP, not the fixed-Tw fold Jacobian.")
end

println("Wrote results to: ", OUTDIR)
end

main()

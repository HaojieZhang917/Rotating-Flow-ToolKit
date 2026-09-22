#!/usr/bin/env julia

using LinearAlgebra
using Random
using Printf
using DelimitedFiles

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKRadialResponse.jl"))
using .BEKIsothermal, .BEKConsistent, .BEKRadialResponse

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT = joinpath(ROOT, "work", "results",
                        "nonsimilar_heated_bek_20260909")

function nonsimilar_mode_residual(state, base, lambda, solution)
    op = solution.operators
    n = length(op.x)
    H = view(state, 1:n)
    F = view(state, n+1:2n)
    G = view(state, 2n+1:3n)
    T = view(state, 3n+1:4n)
    H0 = view(base, 1:n)
    F0 = view(base, n+1:2n)
    G0 = view(base, 2n+1:3n)
    T0 = view(base, 3n+1:4n)
    D, D2 = op.deta, op.deta2
    Feta, Geta, Teta = D * F, D * G, D * T
    xFx = lambda .* (F .- F0)
    xGx = lambda .* (G .- G0)
    xTx = lambda .* (T .- T0)
    Ro, gamma, Pr = solution.Ro, solution.gamma, solution.Pr
    Co = 2 - Ro - Ro^2
    s = Co / 2
    omega = s + Ro
    chi = 1 .- gamma .* (T .- 1)
    q = s .+ Ro .* G

    residual = vcat(
        D * H + 2F + xFx,
        D2 * F + chi .* (Ro .* (F.^2 + F .* xFx + H .* Feta) .-
                          q.^2 ./ Ro) .+ omega^2 / Ro,
        D2 * G + chi .* (Ro .* (2 .* F .* G + F .* xGx +
                                 H .* Geta) .+ Co .* F),
        D2 * T + Pr .* Ro .* (F .* xTx + H .* Teta),
    )

    rH, rF = 1:n, n+1:2n
    rG, rT = 2n+1:3n, 3n+1:4n
    for (row, value) in (
        (first(rH), H[1]),
        (first(rF), F[1]),
        (first(rG), G[1]),
        (last(rF), F[end]),
        (last(rG), G[end] - 1),
        (first(rT), T[1] - solution.Tw),
        (last(rT), T[end] - 1),
    )
        residual[row] = value
    end
    residual[last(rH)] = dot(op.dx[end, :], H)
    residual
end

function main()
    mkpath(OUTPUT)
    rng = MersenneTwister(20260909)
    rows = Vector{Vector{Float64}}()
    for Ro in (-0.50, -0.4703674, -0.40)
        solution = BEKConsistent.solve_consistent_isothermal(
            Ro; degree=60, gamma=1.0, Pr=0.72)
        base = BEKConsistent.pack(solution.fields)
        J, Mr = radial_response_pencil(solution)
        for lambda in (-1.3, -0.2, 0.0, 0.4, 2.0)
            direction = randn(rng, length(base))
            direction ./= norm(direction)
            step = 1e-5
            plus = nonsimilar_mode_residual(
                base .+ step .* direction, base, lambda, solution)
            minus = nonsimilar_mode_residual(
                base .- step .* direction, base, lambda, solution)
            finite_difference = (plus .- minus) ./ (2step)
            analytic = (J .+ lambda .* Mr) * direction
            absolute_error = norm(finite_difference - analytic, Inf)
            relative_error = absolute_error / max(norm(analytic, Inf), eps())
            push!(rows, [Ro, lambda, absolute_error, relative_error])
            @printf("Ro=%+.7f lambda=%+.2f abs=%.3e rel=%.3e\n",
                    Ro, lambda, absolute_error, relative_error)
        end
    end
    open(joinpath(OUTPUT, "radial_pencil_directional_check.csv"), "w") do io
        println(io, "Ro,lambda,absolute_error,relative_error")
        writedlm(io, rows, ',')
    end
end

main()

#!/usr/bin/env julia

"""
Algebraic regression for the non-similar heated BEK boundary-layer model.

This is not a PDE solver. It verifies at randomly generated local states
that zero radial derivatives recover the four continuous similarity
residuals exactly.
"""

using Random
using Printf

function residuals(H, Heta, F, Fx, Feta, Fetaeta,
                   G, Gx, Geta, Getaeta,
                   T, Tx, Teta, Tetaeta, x, Ro, gamma, Pr)
    Co = 2 - Ro - Ro^2
    s = Co / 2
    omega = s + Ro
    chi = 1 - gamma * (T - 1)
    q = s + Ro * G
    Rpde = (
        Heta + 2F + x * Fx,
        Fetaeta + chi * (Ro * (F^2 + x * F * Fx + H * Feta) -
                          q^2 / Ro) + omega^2 / Ro,
        Getaeta + chi * (Ro * (2F * G + x * F * Gx + H * Geta) +
                           Co * F),
        Tetaeta + Pr * Ro * (x * F * Tx + H * Teta),
    )
    Rode = (
        Heta + 2F,
        Fetaeta + chi * (Ro * (F^2 + H * Feta) - q^2 / Ro) +
                          omega^2 / Ro,
        Getaeta + chi * (Ro * (2F * G + H * Geta) + Co * F),
        Tetaeta + Pr * Ro * H * Teta,
    )
    Rpde, Rode
end

function main()
    rng = MersenneTwister(20260909)
    maxerr = 0.0
    for Ro in (-0.50, -0.4703674, -0.40), _ in 1:1000
        H, Heta, F, Feta, Fetaeta, G, Geta, Getaeta,
        T, Teta, Tetaeta = randn(rng, 11)
        x = rand(rng) + 0.05
        Rpde, Rode = residuals(H, Heta, F, 0.0, Feta, Fetaeta,
                               G, 0.0, Geta, Getaeta,
                               T, 0.0, Teta, Tetaeta,
                               x, Ro, 1.0, 0.72)
        maxerr = max(maxerr, maximum(abs.(Rpde .- Rode)))
    end

    @printf("similarity-gate max algebraic residual difference = %.3e\n",
            maxerr)
    maxerr <= 32eps(Float64) ||
        error("non-similar equations do not recover similarity residuals")
end

main()

#!/usr/bin/env julia

using Printf

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const INPUT = joinpath(ROOT,"work","results","corrected_manual_tangent_convergence_20260907")
const DEGREES = [200,240,280,320]
const EPSILONS = [0.02,0.015,0.01]

function rows(path)
    lines = filter(!isempty,strip.(readlines(path)))
    names = Symbol.(split(lines[1],','))
    [NamedTuple{Tuple(names)}(Tuple(parse.(Float64,split(line,',')))) for line in lines[2:end]]
end

allrows = reduce(vcat,[rows(joinpath(INPUT,"N$(n)","tangent_table.csv")) for n in DEGREES])
open(joinpath(INPUT,"convergence_summary.csv"),"w") do io
    println(io,"epsilon,dRo_N200,dRo_N240,dRo_N280,dRo_N320,dTw_N200,dTw_N240,dTw_N280,dTw_N320,relative_dRo_N280_N320,relative_dTw_N280_N320,slope_N320,relative_slope_to_limit")
    for epsilon in EPSILONS
        r = [only(filter(q->q.N==n && isapprox(q.epsilon,epsilon;atol=1e-14),allrows)) for n in DEGREES]
        rdro = abs(r[4].dRo/r[3].dRo-1)
        rdtw = abs(r[4].dTw/r[3].dTw-1)
        slope_limit = -1.021771560631
        @printf(io,"%.12e",epsilon)
        for q in r; @printf(io,",%.12e",q.dRo); end
        for q in r; @printf(io,",%.12e",q.dTw); end
        @printf(io,",%.12e,%.12e,%.12e,%.12e\n",rdro,rdtw,r[4].dTw_dRo,
                abs(r[4].dTw_dRo/slope_limit-1))
    end
end

open(joinpath(INPUT,"scope.txt"),"w") do io
    println(io,"Corrected-energy finite-epsilon augmented-fold tangent convergence.")
    println(io,"N=280 versus N=320 spread is a two-grid diagnostic, not a rigorous error bound.")
    println(io,"No epsilon-to-zero absolute coefficient is reported because dRo and dTw are not uniformly converged at epsilon=0.01.")
    println(io,"Adjoint and direct slope agreement checks the same discrete tangent equation; distance to the analytic limit is reported separately.")
end

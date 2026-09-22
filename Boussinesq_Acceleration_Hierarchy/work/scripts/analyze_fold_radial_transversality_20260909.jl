#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKRadialResponse.jl"))
using .BEKIsothermal, .BEKConsistent, .BEKRadialResponse
using LinearAlgebra
using DelimitedFiles
using Printf

function option(prefix, default)
    for argument in ARGS
        startswith(argument, prefix) &&
            return split(argument, "="; limit=2)[2]
    end
    default
end

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SOURCE = joinpath(ROOT, "work", "results",
    "corrected_energy_epsilon_fold_chain")
const OUTPUT = joinpath(ROOT, "work", "results",
    "fold_radial_transversality_20260909")
const PR = 0.72
const GAMMA = 1.0
const MAP_A = parse(Float64,option("--map-a=","2.0"))
const MAP_B = parse(Float64,option("--map-b=","0.6"))
const MAP_C = parse(Float64,option("--map-c=","0.5"))
const RUN_OVERRIDE = option("--run=", "")
const OUTPUT_SUFFIX = option("--output-suffix=", "")

const RUNS = Dict(
    160 => "convergence_N160",
    200 => "convergence_N200",
    240 => "convergence_N240",
    280 => "production_N280",
    320 => "production_N320",
)

const DEGREES = parse.(Int, split(option("--degrees=", "160,200,240"), ','))
const EPSILONS = parse.(Float64, split(option("--eps=", "0.1,0.05,0.02"), ','))

function numeric_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    names = Symbol.(split(lines[1], ','))
    [NamedTuple{Tuple(names)}(Tuple(parse.(Float64, split(line, ','))))
        for line in lines[2:end]]
end

epsilon_token(epsilon) = replace(@sprintf("%.4g", epsilon),
    "." => "p", "-" => "m")

function load_solution(degree, epsilon)
    run = isempty(RUN_OVERRIDE) ? RUNS[degree] : RUN_OVERRIDE
    directory = joinpath(SOURCE, run)
    table = numeric_csv(joinpath(directory, "epsilon_fold_table.csv"))
    rows = filter(row -> isapprox(row.epsilon, epsilon; atol=1e-13, rtol=0), table)
    length(rows) == 1 || error("epsilon=$epsilon absent or repeated in $run")
    row = only(rows)
    profile = readdlm(joinpath(directory, "profiles",
        "profile_epsilon_$(epsilon_token(epsilon)).csv"), ','; skipstart=1)
    op = BEKIsothermal._operators(degree, MAP_A, MAP_B, MAP_C)
    finite = isfinite.(op.eta)
    maximum(abs.(profile[finite, 1] .- op.eta[finite])) < 1e-7 ||
        error("profile/operator grid mismatch for $run epsilon=$epsilon")
    fields = permutedims(profile[:, 2:5])
    solution = BEKConsistent.ConsistentSolution(row.Ro_f,
        2-row.Ro_f-row.Ro_f^2, GAMMA, PR, op, fields, row.Tw_f,
        row.Hinf, row.residual, 0)
    solution, row
end

function null_pair(J, n)
    rownorm = max.([norm(view(J, row, :)) for row in axes(J, 1)], 1e-14)
    scaled = J ./ rownorm
    factor = svd(scaled)
    v = factor.V[:, end]
    # A mesh-independent amplitude convention for c2: max(abs(v_G)) = 1.
    g = view(v, 2n+1:3n)
    pivot = argmax(abs.(g)) + 2n
    v[pivot] < 0 && (v .*= -1)
    v ./= maximum(abs, g)

    # If u'*(D*J)=0 with D=diag(1/rownorm), then
    # w'=u'*D is the algebraic left null covector of J.
    wraw = factor.U[:, end] ./ rownorm
    wraw ./= norm(wraw)
    pairing = dot(wraw, v)
    abs(pairing) > 1e-12 || error("ill-conditioned left/right pairing")
    w = wraw ./ pairing
    v, w, wraw, pairing, factor.S
end

function weighted_representation_check(w, v, Mr)
    # Any positive diagonal quadrature representation W gives the continuous
    # adjoint samples w_c=W^{-1}w and the identical discrete functional
    # w_c' W Mr v = w' Mr v.  A deliberately nonuniform W checks the bookkeeping.
    m = length(w)
    weights = 0.2 .+ collect(1:m) ./ m
    wc = w ./ weights
    algebraic = dot(w, Mr*v)
    represented = dot(wc, weights .* (Mr*v))
    abs(represented-algebraic) / max(abs(algebraic), 1e-30)
end

function analyze(degree, epsilon)
    solution, source_row = load_solution(degree, epsilon)
    state = BEKConsistent.pack(solution.fields)
    residual, J = BEKConsistent._residual_jacobian(state, solution.Ro,
        solution.Tw, solution.gamma, solution.Pr, solution.operators)
    _, Mr = BEKRadialResponse.radial_response_pencil(solution)
    Drr, Drgeo = BEKRadialResponse.radial_diffusion_operators(solution)
    n = length(solution.operators.x)
    v, w, wraw, raw_wv, singular = null_pair(J, n)

    K = BEKConsistent._jacobian_directional_derivative(state, solution.Ro,
        solution.Tw, solution.gamma, solution.Pr, solution.operators, v)
    rtw = BEKConsistent._tw_derivative(state, solution.operators)
    rro = BEKConsistent._residual_ro_derivative(state, solution.Ro,
        solution.Tw, solution.gamma, solution.Pr, solution.operators)
    Mrv = Mr*v
    raw_cr = dot(wraw, Mrv)
    cr_angle = raw_cr / max(norm(Mrv), 1e-30)
    cr = dot(w, Mrv)
    c2 = 0.5dot(w, K*v)
    ctw = dot(w, rtw)
    cro = dot(w, rro)
    crr = dot(w, Drr*v)
    crgeo = dot(w, Drgeo*v)
    rF=n+1:2n; rG=2n+1:3n; rT=3n+1:4n
    crr_F=dot(view(w,rF),view(v,rF))
    crr_G=dot(view(w,rG),view(v,rG))
    crr_T=dot(view(w,rT),view(v,rT))
    (; N=degree, epsilon, Ro=solution.Ro, Tw=solution.Tw,
       source_residual=source_row.residual,
       residual=norm(residual, Inf),
       sigma_min=singular[end], sigma_ratio=singular[end]/singular[1],
       right_residual=norm(J*v, Inf)/max(norm(v, Inf), 1e-30),
       left_residual=norm(transpose(J)*w, Inf)/max(norm(w, Inf), 1e-30),
       wv=dot(w, v), raw_wv, raw_cr, cr_angle,
       cr, c2, c2_by_cr=c2/cr, ctw, cro,
       cr_by_ctw=cr/ctw, c2_by_ctw=c2/ctw, cro_by_ctw=cro/ctw,
       crr, crgeo, crr_by_ctw=crr/ctw,
       crgeo_by_ctw=crgeo/ctw,
       crr_F_by_ctw=crr_F/ctw, crr_G_by_ctw=crr_G/ctw,
       crr_T_by_ctw=crr_T/ctw,
       weighted_pairing_relative_error=weighted_representation_check(w,v,Mr))
end

function main()
    mkpath(OUTPUT)
    rows = NamedTuple[]
    for degree in DEGREES, epsilon in EPSILONS
        @printf("N=%d epsilon=%.5g ...\n", degree, epsilon)
        row = analyze(degree, epsilon)
        push!(rows, row)
        @printf("  Ro=%+.10f cr=%+.8e c2=%+.8e c2/cr=%+.8e rR=%.2e rL=%.2e\n",
            row.Ro, row.cr, row.c2, row.c2_by_cr,
            row.right_residual, row.left_residual)
    end
    suffix = isempty(OUTPUT_SUFFIX) ? "" : "_$(OUTPUT_SUFFIX)"
    path = joinpath(OUTPUT, "fold_radial_transversality$(suffix).csv")
    open(path, "w") do io
        names = propertynames(first(rows))
        println(io, join(string.(names), ','))
        for row in rows
            for (index, name) in enumerate(names)
                index > 1 && print(io, ',')
                value = getproperty(row, name)
                value isa Integer ? print(io, value) : @printf(io, "%.16e", value)
            end
            println(io)
        end
    end
    open(joinpath(OUTPUT, "scope.txt"), "w") do io
        println(io, "Source: saved corrected-energy finite-epsilon fold profiles.")
        println(io, "No fold continuation or profile correction is performed here.")
        println(io, "Right mode convention: max(abs(vG))=1 and its maximum component is positive.")
        println(io, "Left covector convention: w^T v=1.")
        println(io, "cr is invariant under reciprocal rescaling of the null pair; c2 follows the stated v convention.")
        println(io, "This is a steady axisymmetric radial-response diagnostic, not temporal or 3D stability.")
    end
    fit_rows = filter(row -> row.N >= maximum(DEGREES)-40 &&
               row.epsilon >= 0.01, rows)
    open(joinpath(OUTPUT, "scaling_fits$(suffix).csv"), "w") do io
        println(io, "N,quantity,prefactor,power,points")
        for degree in sort(unique(getfield.(fit_rows, :N)))
            subset = sort(filter(row -> row.N == degree, fit_rows);
                by=row -> row.epsilon)
            for name in (:cr_by_ctw, :c2_by_ctw, :c2_by_cr,
                         :crr_by_ctw, :crgeo_by_ctw)
                usable = filter(row -> isfinite(getproperty(row,name)) &&
                    getproperty(row,name) != 0, subset)
                x = log.(getfield.(usable,:epsilon))
                y = log.(abs.(getproperty.(usable,name)))
                coefficient = hcat(ones(length(x)),x) \ y
                @printf(io,"%d,%s,%.12e,%.12e,%d\n",degree,string(name),
                    exp(coefficient[1]),coefficient[2],length(x))
            end
        end
    end
    println("Output: $path")
end

main()

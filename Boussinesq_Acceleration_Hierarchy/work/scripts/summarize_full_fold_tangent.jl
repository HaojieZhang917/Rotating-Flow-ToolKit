#!/usr/bin/env julia

using Printf

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const INPUT = joinpath(ROOT,"work","results","full_fold_tangent_analysis")
const DEGREES = [240,280,320]

function numeric_csv(path)
    lines = filter(!isempty,strip.(readlines(path)))
    header = Symbol.(split(lines[1],','))
    [NamedTuple{Tuple(header)}(Tuple(parse.(Float64,split(line,',')))) for line in lines[2:end]]
end

function main()
    rows = NamedTuple[]
    for degree in DEGREES
        append!(rows,numeric_csv(joinpath(INPUT,"N$(degree)","tangent_table.csv")))
    end
    sort!(rows,by=r->(r.epsilon,r.N),rev=true)

    open(joinpath(INPUT,"resolution_comparison.csv"),"w") do io
        println(io,"N,epsilon,dRo_depsilon,dTw_depsilon,dTw_dRo,left_slope,tangent_residual,fold_quadratic,Ro,Tw")
        for r in rows
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                Int(r.N),r.epsilon,r.dRo,r.dTw,r.dTw_dRo,r.left_slope,
                r.tangent_residual,r.fold_quadratic,r.Ro,r.Tw)
        end
    end

    open(joinpath(INPUT,"right_mode_outer_scaling.csv"),"w") do io
        println(io,"N,epsilon,Hmax_over_Gmax_epsilon,Fmax_over_Gmax_epsilon2,Tmax_over_Gmax,abs_null_Hinf_over_Gmax_epsilon")
        for r in rows
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e\n",Int(r.N),r.epsilon,
                r.Hmax/(r.Gmax*r.epsilon),r.Fmax/(r.Gmax*r.epsilon^2),
                r.Tmax/r.Gmax,abs(r.null_Hinf)/(r.Gmax*r.epsilon))
        end
    end

    extrapolations = NamedTuple[]
    for degree in DEGREES
        r20 = only(filter(r->Int(r.N)==degree && r.epsilon==0.02,rows))
        r10 = only(filter(r->Int(r.N)==degree && r.epsilon==0.01,rows))
        # If d(parameter)/d(epsilon)=A+2B*epsilon+..., these two points
        # give the first-order epsilon->0 intercept without using the
        # resolution-contaminated epsilon=0.005 tangent.
        AR = 2r10.dRo-r20.dRo
        AT = 2r10.dTw-r20.dTw
        push!(extrapolations,(;N=degree,AR,AT,ratio=AT/AR))
    end
    open(joinpath(INPUT,"resolved_tangent_limit_estimates.csv"),"w") do io
        println(io,"N,A_R_from_epsilon_0p02_0p01,A_T_from_epsilon_0p02_0p01,A_T_over_A_R")
        for r in extrapolations
            @printf(io,"%d,%.12e,%.12e,%.12e\n",r.N,r.AR,r.AT,r.ratio)
        end
    end

    best = only(filter(r->r.N==320,extrapolations))
    fit_AR=-0.420342495108; fit_AT=0.510660114111
    open(joinpath(INPUT,"summary.txt"),"w") do io
        println(io,"method=implicit derivative of complete finite-epsilon augmented fold system")
        println(io,"resolved_limit_estimate=linear extrapolation of tangent at epsilon 0.02 and 0.01")
        println(io,"epsilon_0p005_status=resolution contaminated; excluded from limit estimate")
        @printf(io,"N320_A_R=%.12f\n",best.AR)
        @printf(io,"N320_A_T=%.12f\n",best.AT)
        @printf(io,"N320_ratio=%.12f\n",best.ratio)
        @printf(io,"fit_A_R=%.12f\n",fit_AR)
        @printf(io,"fit_A_T=%.12f\n",fit_AT)
        @printf(io,"A_R_difference=%.12e\n",best.AR-fit_AR)
        @printf(io,"A_T_difference=%.12e\n",best.AT-fit_AT)
        @printf(io,"A_R_relative_difference=%.12e\n",(best.AR-fit_AR)/fit_AR)
        @printf(io,"A_T_relative_difference=%.12e\n",(best.AT-fit_AT)/fit_AT)
    end
    @printf("N320 resolved tangent intercept: A_R=%+.9f A_T=%+.9f ratio=%+.9f\n",
            best.AR,best.AT,best.ratio)
end

main()

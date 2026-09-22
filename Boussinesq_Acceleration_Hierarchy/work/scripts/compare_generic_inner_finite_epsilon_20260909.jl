#!/usr/bin/env julia

using LinearAlgebra
using Printf

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const FULL=joinpath(ROOT,"work","results","generic_tail_zero_core_20260909",
    "N320_a2p0000_b0p6000_c0p5000")
const INNER=joinpath(ROOT,"work","results","generic_inner_first_order_20260909")
const OUT=joinpath(INNER,"finite_epsilon_comparison.csv")

function table(path)
    lines=filter(!isempty,strip.(readlines(path))); h=split(first(lines),',')
    [Dict(h .=> split(line,',')) for line in lines[2:end]]
end
num(r,k)=parse(Float64,r[k])
token(x)=replace(@sprintf("%.4f",abs(x)),"."=>"p")
inner_token(x)=replace(@sprintf("%.6f",abs(x)),"."=>"p")
sfun(r)=(2-r-r^2)/2
tcurve(r)=1+(1-((sfun(r)+r)/sfun(r))^2)

function interpolate(x,y,grid)
    out=similar(grid)
    j=1
    for (k,z) in pairs(grid)
        while j<length(x)-1 && x[j+1]<z; j+=1; end
        t=(z-x[j])/(x[j+1]-x[j]); out[k]=(1-t)*y[j]+t*y[j+1]
    end
    out
end

function curve(rows,xkey,ykey,grid)
    x=[num(r,xkey) for r in rows]; y=[num(r,ykey) for r in rows]
    interpolate(x,y,grid)
end

function main()
    grid=collect(0.0:.01:8.0)
    ros=[-.469,-.467,-.465,-.46,-.44,-.42,-.40,-.35,-.30]
    open(OUT,"w") do io
        println(io,"Ro,epsilon,field,relative_L2,absolute_rms,max_abs,eta_range")
        for ro in ros
            inner=table(joinpath(INNER,"profile_Ro_$(inner_token(ro)).csv"))
            asym=Dict(
                "F1"=>curve(inner,"eta","F1",grid),
                "H1"=>curve(inner,"eta","H1",grid),
                "G1"=>curve(inner,"eta","G1",grid),
                "Theta1"=>curve(inner,"eta","Theta1",grid))
            for epsilon in (.03,.02,.015,.01)
                full=table(joinpath(FULL,"profiles","Ro_m$(token(ro))",
                    "epsilon_$(token(epsilon)).csv"))
                scaled=Dict(
                    "F1"=>curve(full,"eta","F",grid)/epsilon,
                    "H1"=>curve(full,"eta","H",grid)/epsilon,
                    "G1"=>curve(full,"eta","G",grid)/epsilon,
                    "Theta1"=>(curve(full,"eta","T",grid).-tcurve(ro))/epsilon)
                for field in ("F1","H1","G1","Theta1")
                    error=scaled[field]-asym[field]
                    rel=norm(error)/max(norm(asym[field]),eps(Float64))
                    @printf(io,"%.12e,%.12e,%s,%.12e,%.12e,%.12e,0:8\n",
                        ro,epsilon,field,rel,norm(error)/sqrt(length(error)),maximum(abs,error))
                end
            end
        end
    end
    println("Output: $OUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

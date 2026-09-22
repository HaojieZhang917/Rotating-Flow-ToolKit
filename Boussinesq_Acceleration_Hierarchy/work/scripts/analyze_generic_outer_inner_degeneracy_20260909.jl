#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Plots

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const INNER=joinpath(ROOT,"work","results","generic_inner_first_order_20260909")
const OUTER=joinpath(ROOT,"work","results","generic_outer_family_20260909")
const OUT=joinpath(ROOT,"work","results","generic_outer_inner_degeneracy_20260909")
const RT=-0.470367416464038

function rows(path)
    lines=filter(!isempty,strip.(readlines(path))); h=split(first(lines),',')
    [Dict(h .=> split(line,',')) for line in lines[2:end]]
end
num(r,k)=parse(Float64,r[k])
tok(r)=replace(@sprintf("%.6f",abs(r)),"."=>"p")

function mmer(path; etamax=8.0)
    data=rows(path); eta=[num(r,"eta") for r in data]
    F=[num(r,"F1") for r in data]; H=[num(r,"H1") for r in data]
    keep=findall(x->x<=etamax+1e-12,eta)
    e=eta[keep]; q=(F[keep].^2+H[keep].^2)
    sqrt(sum((q[1:end-1]+q[2:end]).*(e[2:end]-e[1:end-1])/2))
end

function main()
    mkpath(OUT)
    iconv=filter(r->abs(num(r,"L")-24)<1e-10 && abs(num(r,"deta")-.005)<1e-10,
        rows(joinpath(INNER,"convergence.csv")))
    orows=filter(r->round(Int,num(r,"g_steps"))==16000,
        rows(joinpath(OUTER,"family_convergence.csv")))
    odict=Dict(round(num(r,"Ro"),digits=9)=>r for r in orows)
    data=NamedTuple[]
    for r in iconv
        ro=num(r,"Ro"); abs(ro-RT)<=.00301 || continue
        outer=odict[round(ro,digits=9)]
        mm=mmer(joinpath(INNER,"profile_Ro_$(tok(ro)).csv"))
        push!(data,(;ro,mu=ro-RT,h0=num(outer,"h0"),D1=num(r,"D1"),
            tau1=num(r,"tau1"),Mmer1=mm))
    end
    sort!(data,by=x->x.ro)
    open(joinpath(OUT,"local_degeneracy_data.csv"),"w") do io
        println(io,"Ro,mu,h0,D1,tau1,Mmer1_eta0to8")
        for x in data
            @printf(io,"%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                x.ro,x.mu,x.h0,x.D1,x.tau1,x.Mmer1)
        end
    end
    open(joinpath(OUT,"local_fit_summary.csv"),"w") do io
        println(io,"quantity,model,leading_coefficient,rms_residual,npoints")
        for quantity in (:h0,:D1,:tau1)
            y=getfield.(data,quantity); X=hcat(getfield.(data,:mu),getfield.(data,:mu).^2)
            c=X\y; rms=norm(y-X*c)/sqrt(length(y))
            @printf(io,"%s,c1*mu+c2*mu2,%.12e,%.12e,%d\n",String(quantity),c[1],rms,length(y))
        end
        y=getfield.(data,:Mmer1); x=abs.(getfield.(data,:mu)); X=hcat(x,x.^2)
        c=X\y; rms=norm(y-X*c)/sqrt(length(y))
        @printf(io,"Mmer1,c1*abs(mu)+c2*mu2,%.12e,%.12e,%d\n",c[1],rms,length(y))
    end
    p=plot(;xlabel="Ro-Ro_t",ylabel="scaled coefficient",framestyle=:box,
        gridalpha=.25,size=(820,500))
    plot!(getfield.(data,:mu),getfield.(data,:D1),marker=:circle,label="D1")
    plot!(getfield.(data,:mu),-getfield.(data,:h0),marker=:square,label="-h0")
    plot!(getfield.(data,:mu),getfield.(data,:Mmer1),marker=:diamond,label="Mmer1")
    vline!([0.0],linestyle=:dash,color=:black,label="Ro_t")
    savefig(joinpath(OUT,"local_linear_degeneracy.png"))
    println("Output: $OUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

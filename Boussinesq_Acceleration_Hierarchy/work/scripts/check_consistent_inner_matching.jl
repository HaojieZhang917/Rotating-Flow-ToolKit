#!/usr/bin/env julia
# Existing-data diagnostics only: no solver, continuation, or fitted parameters.
using DelimitedFiles, Printf

const ROOT = normpath(joinpath(@__DIR__,"..",".."))
const OUT = joinpath(ROOT,"work","results","consistent_inner_leading_matching")
const DATA = joinpath(ROOT,"work","results","consistent_epsilon_fold_chain",
                      "convergence_N320")

function algebra(r,tw;gamma=1.0)
    s=(2-r-r^2)/2
    omega=s+r
    chi=1-gamma*(tw-1)
    chi>0 || error("positive chi required")
    gm=(omega/sqrt(chi)-s)/r
    tw_zero=1+(1-(omega/s)^2)/gamma
    (s=s,omega=omega,chi=chi,gm=gm,tw_zero=tw_zero,
     defect=tw-tw_zero,
     balance=chi*(s+r*gm)^2-omega^2)
end

function at_eta(profile,eta)
    x=profile[:,1]
    k=searchsortedlast(x,eta)
    1<=k<length(x) || error("eta outside finite sample interval")
    weight=(eta-x[k])/(x[k+1]-x[k])
    (1-weight).*profile[k,2:5]+weight.*profile[k+1,2:5]
end

function main()
    mkpath(OUT)
    fits=split.(readlines(joinpath(ROOT,"work","results",
        "consistent_fold_tail_scaling","polynomial_fits.csv"))[2:end],',')
    open(joinpath(OUT,"limit_compatibility.csv"),"w") do io
        println(io,"window,Ro_t,Tw_t,g_m,Tw_zero_core,Tw_defect,radial_balance")
        for window in ("all_five","drop_epsilon_0p1")
            r=parse(Float64,only([a for a in fits if a[1]=="Ro" && a[2]==window && a[3]=="2"])[4])
            tw=parse(Float64,only([a for a in fits if a[1]=="Tw" && a[2]==window && a[3]=="2"])[4])
            d=algebra(r,tw)
            @printf(io,"%s,%.12f,%.12f,%.12e,%.12f,%.12e,%.12e\n",
                window,r,tw,d.gm,d.tw_zero,d.defect,d.balance)
            @printf("%s: gm=%+.7e; Tw - Tw_zero_core=%+.7e\n",window,d.gm,d.defect)
            abs(d.balance)<1e-12 || error("radial balance check failed")
        end
    end

    rows=split.(readlines(joinpath(DATA,"epsilon_fold_table.csv"))[2:end],',')
    open(joinpath(OUT,"fixed_eta_profiles.csv"),"w") do io
        println(io,"epsilon,eta,H,F,G,T,T_minus_Tw,H_over_epsilon,G_over_epsilon,Tdrop_over_epsilon")
        for row in rows
            eps=parse(Float64,row[1]); tw=parse(Float64,row[3])
            token=replace(@sprintf("%.4g",eps),"."=>"p","-"=>"m")
            profile=readdlm(joinpath(DATA,"profiles","profile_epsilon_$(token).csv"),',';skipstart=1)
            for eta in (1.0,2.0,5.0)
                H,F,G,T=at_eta(profile,eta)
                @printf(io,"%.6g,%.1f,%.12e,%.12e,%.12e,%.12f,%.12e,%.12e,%.12e,%.12e\n",
                    eps,eta,H,F,G,T,T-tw,H/eps,G/eps,(tw-T)/eps)
                eta==1 && @printf("eps=%.3g eta=1: H=%+.4e F=%+.4e G=%+.4e Tw-T=%.4e\n",
                    eps,H,F,G,tw-T)
            end
        end
    end
    open(joinpath(OUT,"scope.txt"),"w") do io
        println(io,"No new boundary-value problem or continuation was solved.")
        println(io,"Fields at eta=1,2,5 use piecewise linear interpolation of saved N=320 profiles.")
        println(io,"Fixed-eta samples are not overlap plateau measurements.")
        println(io,"g_m values use extrapolated limiting coordinates and inherit their uncertainty.")
        println(io,"Near-zero g_m supports, but does not prove, a trivial leading inner velocity field.")
    end
end

main()

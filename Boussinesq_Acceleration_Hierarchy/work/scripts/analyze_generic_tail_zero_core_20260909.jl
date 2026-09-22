#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Plots

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const INPUT = joinpath(ROOT, "work", "results", "generic_tail_zero_core_20260909")
const OUTPUT = joinpath(INPUT, "summary")

function rows(path)
    lines = filter(!isempty, strip.(readlines(path)))
    header = split(first(lines), ',')
    [Dict(header .=> split(line, ',')) for line in lines[2:end]]
end
num(row, key) = parse(Float64, row[key])

function fit_intercept(epsilon, values, order)
    X = hcat((epsilon.^k for k in 0:order)...)
    coefficient = X \ values
    residual = norm(values-X*coefficient)/sqrt(length(values))
    coefficient[1], residual
end

function root_linear(x, y)
    c = hcat(ones(length(x)), x) \ y
    -c[1]/c[2], c
end

function main()
    mkpath(OUTPUT)
    dirs = filter(d -> isfile(joinpath(INPUT, d, "diagnostics.csv")), readdir(INPUT))
    coefficient_rows = NamedTuple[]
    for dir in sort(dirs)
        data = rows(joinpath(INPUT, dir, "diagnostics.csv"))
        haskey(first(data), "Mmer_over_epsilon_eta0to8") || continue
        degree = round(Int, num(first(data), "N"))
        for ro in sort(unique(num(r,"Ro") for r in data))
            group = sort(filter(r -> num(r,"Ro") == ro, data), by=r -> num(r,"epsilon"))
            for (window, cutoff, order) in (("eps_le_0p03_quad",0.0301,2),
                                             ("eps_le_0p02_quad",0.0201,2),
                                             ("eps_le_0p02_linear",0.0201,1))
                selected = filter(r -> num(r,"epsilon") <= cutoff, group)
                length(selected) >= order+1 || continue
                eps = [num(r,"epsilon") for r in selected]
                quantities = (
                    D1=[num(r,"D_algebra")/num(r,"epsilon") for r in selected],
                    h0=[num(r,"H_over_epsilon_eta8") for r in selected],
                    h0_integral=[num(r,"minus2intF_over_epsilon_eta0to8") for r in selected],
                    Mmer1=[num(r,"Mmer_over_epsilon_eta0to8") for r in selected])
                for (name, values) in pairs(quantities)
                    intercept, residual = fit_intercept(eps, values, order)
                    push!(coefficient_rows, (; run=dir, degree, ro, window,
                        quantity=String(name), intercept, residual,
                        npoints=length(selected)))
                end
            end
            if length(group) >= 2
                a, b = group[1], group[2]
                e1, e2 = num(a,"epsilon"), num(b,"epsilon")
                m1 = e1*num(a,"Mmer_over_epsilon_eta0to8")
                m2 = e2*num(b,"Mmer_over_epsilon_eta0to8")
                p = log(abs(m2/m1))/log(e2/e1)
                push!(coefficient_rows, (; run=dir, degree, ro, window="smallest_pair",
                    quantity="p_Mmer", intercept=p, residual=NaN, npoints=2))
            end
        end
    end

    coeffile = joinpath(OUTPUT, "generic_tail_first_order_coefficients.csv")
    open(coeffile,"w") do io
        println(io,"run,N,Ro,window,quantity,intercept,rms_residual,npoints")
        for r in coefficient_rows
            @printf(io,"%s,%d,%.12e,%s,%s,%.12e,%.12e,%d\n",
                r.run,r.degree,r.ro,r.window,r.quantity,r.intercept,r.residual,r.npoints)
        end
    end

    roots = NamedTuple[]
    for run in sort(unique(r.run for r in coefficient_rows)),
        degree in sort(unique(r.degree for r in coefficient_rows if r.run==run)),
        window in ("eps_le_0p03_quad","eps_le_0p02_quad","eps_le_0p02_linear"),
        quantity in ("D1","h0","h0_integral")
        selected = filter(r -> r.run==run && r.degree==degree && r.window==window &&
            r.quantity==quantity && -0.4703 <= r.ro <= -0.4589, coefficient_rows)
        length(selected) >= 2 || continue
        sort!(selected, by=r->r.ro)
        root, c = root_linear([r.ro for r in selected], [r.intercept for r in selected])
        push!(roots,(;run,degree,window,quantity,root,slope=c[2],npoints=length(selected)))
    end
    open(joinpath(OUTPUT,"critical_ro_from_first_order_zeros.csv"),"w") do io
        println(io,"run,N,window,quantity,Ro_zero,slope,npoints")
        for r in roots
            @printf(io,"%s,%d,%s,%s,%.12e,%.12e,%d\n",
                r.run,r.degree,r.window,r.quantity,r.root,r.slope,r.npoints)
        end
    end

    fixed_coefficients = NamedTuple[]
    for dir in sort(dirs)
        path = joinpath(INPUT,dir,"fixed_eta.csv")
        isfile(path) || continue
        data = rows(path)
        degree = round(Int,num(first(data),"N"))
        for ro in sort(unique(num(r,"Ro") for r in data)),
            eta in sort(unique(num(r,"eta") for r in data)),
            (window,cutoff,order) in (("eps_le_0p03_quad",.0301,2),
                                      ("eps_le_0p02_quad",.0201,2))
            selected = sort(filter(r -> num(r,"Ro")==ro && num(r,"eta")==eta &&
                num(r,"epsilon")<=cutoff,data),by=r->num(r,"epsilon"))
            length(selected)>=order+1 || continue
            eps = [num(r,"epsilon") for r in selected]
            for quantity in ("H","F")
                values = [num(r,quantity)/num(r,"epsilon") for r in selected]
                intercept,residual=fit_intercept(eps,values,order)
                push!(fixed_coefficients,(;run=dir,degree,ro,eta,window,quantity,
                    intercept,residual,npoints=length(selected)))
            end
        end
    end
    open(joinpath(OUTPUT,"fixed_eta_first_order_coefficients.csv"),"w") do io
        println(io,"run,N,Ro,eta,window,quantity,intercept,rms_residual,npoints")
        for r in fixed_coefficients
            @printf(io,"%s,%d,%.12e,%.4f,%s,%s,%.12e,%.12e,%d\n",r.run,r.degree,
                r.ro,r.eta,r.window,r.quantity,r.intercept,r.residual,r.npoints)
        end
    end
    open(joinpath(OUTPUT,"critical_ro_from_fixed_eta_H1.csv"),"w") do io
        println(io,"run,N,eta,window,Ro_zero,slope,npoints")
        for run in sort(unique(r.run for r in fixed_coefficients)),
            eta in sort(unique(r.eta for r in fixed_coefficients if r.run==run)),
            window in ("eps_le_0p03_quad","eps_le_0p02_quad")
            selected=filter(r->r.run==run && r.eta==eta && r.window==window &&
                r.quantity=="H" && -.4703<=r.ro<=-.4589,fixed_coefficients)
            length(selected)>=2 || continue
            root,c=root_linear([r.ro for r in selected],[r.intercept for r in selected])
            @printf(io,"%s,%d,%.4f,%s,%.12e,%.12e,%d\n",run,first(selected).degree,
                eta,window,root,c[2],length(selected))
        end
    end

    for quantity in ("D1","h0","Mmer1")
        plot(;xlabel="Rossby number",ylabel=quantity,framestyle=:box,
            gridalpha=.25,legend=:best,size=(820,500))
        for run in sort(unique(r.run for r in coefficient_rows))
            selected = filter(r -> r.run==run && r.window=="eps_le_0p03_quad" &&
                r.quantity==quantity, coefficient_rows)
            sort!(selected,by=r->r.ro)
            isempty(selected) || plot!([r.ro for r in selected],
                [r.intercept for r in selected],marker=:circle,label=run)
        end
        vline!([-0.4703674163],linestyle=:dash,color=:black,label="outer Ro_t")
        savefig(joinpath(OUTPUT,"$(quantity)_versus_Ro.png"))
    end
    println("Output: $OUTPUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

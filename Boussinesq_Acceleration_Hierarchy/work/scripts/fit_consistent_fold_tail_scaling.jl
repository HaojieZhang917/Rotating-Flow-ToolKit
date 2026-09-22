#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Plots

const INPUT = joinpath(@__DIR__,"..","results",
    "consistent_epsilon_fold_chain","recommended_epsilon_fold_table.csv")
const OUTPUT = joinpath(@__DIR__,"..","results",
                        "consistent_fold_tail_scaling")

function read_input(path)
    lines = readlines(path)
    rows = NamedTuple[]
    for line in lines[2:end]
        values = split(line,',')
        push!(rows,(epsilon=parse(Float64,values[1]),
            Ro=parse(Float64,values[2]),Tw=parse(Float64,values[3]),
            Ro_error=parse(Float64,values[4]),
            Tw_error=parse(Float64,values[5])))
    end
    sort(rows;by=row->row.epsilon,rev=true)
end

function polynomial_fit(rows,field,order)
    epsilon = getfield.(rows,:epsilon)
    values = getfield.(rows,field)
    design = hcat([epsilon.^power for power in 0:order]...)
    coefficients = design\values
    residual = values-design*coefficients
    (coefficients=coefficients,residual=residual,
     rms=sqrt(sum(abs2,residual)/length(residual)),
     maximum=maximum(abs.(residual)))
end

function effective_exponent(rows,field)
    selected = rows[end-2:end]
    epsilon = getfield.(selected,:epsilon)
    values = getfield.(selected,field)
    ratio = (values[1]-values[2])/(values[2]-values[3])
    exponent = log(abs(ratio))/log(epsilon[1]/epsilon[2])
    (exponent=exponent,ratio=ratio,
     delta_large=values[1]-values[2],
     delta_small=values[2]-values[3])
end

function power_profile_fit(rows,field;pmin=0.5,pmax=1.5,points=20001)
    epsilon = getfield.(rows,:epsilon)
    values = getfield.(rows,field)
    best = nothing
    for p in range(pmin,pmax;length=points)
        design = hcat(ones(length(epsilon)),epsilon.^p)
        coefficients = design\values
        residual = values-design*coefficients
        sse = sum(abs2,residual)
        if best === nothing || sse < best.sse
            best = (p=p,limit=coefficients[1],amplitude=coefficients[2],
                    sse=sse,rms=sqrt(sse/length(values)))
        end
    end
    best
end

function corner_sensitivity(rows,field,error_field;order=2)
    epsilon = getfield.(rows,:epsilon)
    values = getfield.(rows,field)
    errors = getfield.(rows,error_field)
    design = hcat([epsilon.^power for power in 0:order]...)
    coefficients = Vector{Vector{Float64}}()
    for mask in 0:(2^length(rows)-1)
        signs = [((mask >> (index-1)) & 1) == 1 ? 1.0 : -1.0
                 for index in eachindex(rows)]
        push!(coefficients,design\(values.+signs.*errors))
    end
    matrix = hcat(coefficients...)
    (minimum=vec(minimum(matrix;dims=2)),
     maximum=vec(maximum(matrix;dims=2)))
end

function exponent_corner_sensitivity(rows,field,error_field)
    selected = rows[end-2:end]
    values = getfield.(selected,field)
    errors = getfield.(selected,error_field)
    exponents = Float64[]
    for mask in 0:7
        signs = [((mask >> (index-1)) & 1) == 1 ? 1.0 : -1.0
                 for index in 1:3]
        perturbed = values.+signs.*errors
        ratio = (perturbed[1]-perturbed[2])/
                (perturbed[2]-perturbed[3])
        push!(exponents,log(abs(ratio))/log(2))
    end
    (minimum=minimum(exponents),maximum=maximum(exponents))
end

function main()
    mkpath(OUTPUT)
    rows = read_input(INPUT)
    windows = [("all_five",rows),("drop_epsilon_0p1",rows[2:end])]
    fits = NamedTuple[]
    for (label,window) in windows, field in (:Ro,:Tw), order in (1,2)
        fit = polynomial_fit(window,field,order)
        coefficients = vcat(fit.coefficients,
                             fill(NaN,3-length(fit.coefficients)))
        push!(fits,(variable=String(field),window=label,order=order,
            limit=coefficients[1],A=coefficients[2],B=coefficients[3],
            rms=fit.rms,maximum=fit.maximum,n=length(window)))
    end
    open(joinpath(OUTPUT,"polynomial_fits.csv"),"w") do io
        println(io,"variable,window,order,limit,A_linear,B_quadratic,rms_residual,max_residual,npoints")
        for row in fits
            @printf(io,"%s,%s,%d,%.12f,%.12f,%.12f,%.12e,%.12e,%d\n",
                    row.variable,row.window,row.order,row.limit,row.A,row.B,
                    row.rms,row.maximum,row.n)
        end
    end

    exponents = Dict(field=>effective_exponent(rows,field)
                     for field in (:Ro,:Tw))
    profiles = Dict(field=>power_profile_fit(rows[2:end],field)
                    for field in (:Ro,:Tw))
    open(joinpath(OUTPUT,"exponent_diagnostics.csv"),"w") do io
        println(io,"variable,geometric_triple_exponent,difference_ratio,profile_fit_exponent,profile_limit,profile_amplitude,profile_rms")
        for field in (:Ro,:Tw)
            e,p = exponents[field],profiles[field]
            @printf(io,"%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12e\n",
                    field,e.exponent,e.ratio,p.p,p.limit,p.amplitude,p.rms)
        end
    end

    sensitivities = Dict{Symbol,Any}()
    open(joinpath(OUTPUT,"uncertainty_sensitivity.csv"),"w") do io
        println(io,"variable,p_min,p_max,limit_min,limit_max,A_min,A_max,B_min,B_max")
        for (field,error_field) in ((:Ro,:Ro_error),(:Tw,:Tw_error))
            exponent_range = exponent_corner_sensitivity(rows,field,error_field)
            coefficient_range = corner_sensitivity(rows[2:end],field,
                                                     error_field;order=2)
            sensitivities[field] = (p=exponent_range,c=coefficient_range)
            @printf(io,"%s,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f,%.12f\n",
                    field,exponent_range.minimum,exponent_range.maximum,
                    coefficient_range.minimum[1],coefficient_range.maximum[1],
                    coefficient_range.minimum[2],coefficient_range.maximum[2],
                    coefficient_range.minimum[3],coefficient_range.maximum[3])
        end
    end

    # The primary coefficient estimate is the quadratic fit on the four
    # smallest epsilon values; the all-five fit is retained as the requested
    # window-stability check.
    selected = Dict{Symbol,Any}()
    normalized = NamedTuple[]
    for field in (:Ro,:Tw)
        fit = polynomial_fit(rows[2:end],field,2)
        limit,A,B = fit.coefficients
        selected[field] = (limit=limit,A=A,B=B,fit=fit)
        for row in rows
            epsilon = row.epsilon
            value = getfield(row,field)
            leading = (value-limit)/epsilon
            second = (value-limit-A*epsilon)/epsilon^2
            push!(normalized,(variable=String(field),epsilon=epsilon,
                value=value,leading_ratio=leading,
                quadratic_ratio=second))
        end
    end
    open(joinpath(OUTPUT,"normalized_scaling.csv"),"w") do io
        println(io,"variable,epsilon,value,(value-limit)/epsilon,(value-limit-A*epsilon)/epsilon^2")
        for row in normalized
            @printf(io,"%s,%.12e,%.12f,%.12f,%.12f\n",
                    row.variable,row.epsilon,row.value,row.leading_ratio,
                    row.quadratic_ratio)
        end
    end

    open(joinpath(OUTPUT,"scaling_conclusion.txt"),"w") do io
        println(io,"input=existing five-point recommended epsilon-fold table only")
        println(io,"no_new_base_flow_solutions=true")
        for field in (:Ro,:Tw)
            result = selected[field]
            allfit = polynomial_fit(rows,field,2)
            @printf(io,"%s_effective_exponent=%.12f\n",field,
                    exponents[field].exponent)
            @printf(io,"%s_profile_exponent_four_smallest=%.12f\n",field,
                    profiles[field].p)
            @printf(io,"%s_limit_four_smallest_quadratic=%.12f\n",field,
                    result.limit)
            @printf(io,"%s_A_four_smallest_quadratic=%.12f\n",field,
                    result.A)
            @printf(io,"%s_B_four_smallest_quadratic=%.12f\n",field,
                    result.B)
            @printf(io,"%s_limit_all_five_quadratic=%.12f\n",field,
                    allfit.coefficients[1])
            @printf(io,"%s_A_all_five_quadratic=%.12f\n",field,
                    allfit.coefficients[2])
            @printf(io,"%s_effective_exponent_corner_range=[%.12f,%.12f]\n",
                    field,sensitivities[field].p.minimum,
                    sensitivities[field].p.maximum)
            @printf(io,"%s_A_corner_range=[%.12f,%.12f]\n",field,
                    sensitivities[field].c.minimum[2],
                    sensitivities[field].c.maximum[2])
        end
        println(io,"supported_scaling=p_equals_q_equals_1")
        println(io,"scope=numerical_fit_not_matched_asymptotic_derivation")
    end

    default(;framestyle=:box,gridalpha=0.25,linewidth=2,
            guidefontsize=10,tickfontsize=8,legendfontsize=8)
    epsilon_plot = range(0,maximum(getfield.(rows,:epsilon));length=300)
    panels = Any[]
    for field in (:Ro,:Tw)
        result = selected[field]
        curve = result.limit .+ result.A.*epsilon_plot .+
                result.B.*epsilon_plot.^2
        ylabel = field == :Ro ? "Ro_f" : "T_w,f"
        panel = scatter(getfield.(rows,:epsilon),getfield.(rows,field);
            xlabel="epsilon=-Hinf",ylabel=ylabel,label="computed")
        plot!(panel,epsilon_plot,curve;label="quadratic fit, four smallest")
        push!(panels,panel)
    end
    savefig(plot(panels...;layout=(1,2),size=(920,360)),
            joinpath(OUTPUT,"scaling_fits.png"))

    residual_panels = Any[]
    for field in (:Ro,:Tw)
        result = selected[field]
        epsilon = getfield.(rows[2:end],:epsilon)
        values = getfield.(rows[2:end],field)
        first_residual = abs.(values.-result.limit)
        second_residual = abs.(values.-result.limit.-result.A.*epsilon)
        panel = plot(epsilon,first_residual;xscale=:log10,yscale=:log10,
            marker=:circle,label="|y-y_t|",xlabel="epsilon",ylabel="magnitude")
        plot!(panel,epsilon,second_residual;marker=:diamond,
              label="|y-y_t-A epsilon|")
        reference1 = first_residual[end].*(epsilon./epsilon[end])
        reference2 = second_residual[end].*(epsilon./epsilon[end]).^2
        plot!(panel,epsilon,reference1;linestyle=:dash,color=:gray,
              label="O(epsilon)")
        plot!(panel,epsilon,reference2;linestyle=:dot,color=:black,
              label="O(epsilon^2)")
        title!(panel,String(field))
        push!(residual_panels,panel)
    end
    savefig(plot(residual_panels...;layout=(1,2),size=(960,380)),
            joinpath(OUTPUT,"scaling_residuals.png"))

    for field in (:Ro,:Tw)
        result = selected[field]
        @printf("%s: p_eff=%.8f, limit=%.9f, A=%.9f, B=%.9f\n",
                field,exponents[field].exponent,result.limit,result.A,
                result.B)
        @printf("    corner ranges: p=[%.6f,%.6f], A=[%.6f,%.6f]\n",
                sensitivities[field].p.minimum,
                sensitivities[field].p.maximum,
                sensitivities[field].c.minimum[2],
                sensitivities[field].c.maximum[2])
    end
    println("Outputs: $OUTPUT")
end

main()

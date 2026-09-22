#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKRadialResponse.jl"))
using .BEKIsothermal, .BEKConsistent, .BEKRadialResponse
using LinearAlgebra
using DelimitedFiles
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEGREE = 200
const EPSILONS = [0.1, 0.05, 0.02]
const TARGET_AMPLITUDES = [0.015, 0.010, 0.005]
const SOURCE = joinpath(ROOT, "work", "results",
    "corrected_energy_epsilon_fold_chain", "convergence_N200")
const OUTPUT = joinpath(ROOT, "work", "results",
    "fold_radial_transversality_20260909", "normal_form_validation_N200.csv")

function numeric_csv(path)
    lines = filter(!isempty, strip.(readlines(path)))
    names = Symbol.(split(lines[1], ','))
    [NamedTuple{Tuple(names)}(Tuple(parse.(Float64, split(line, ','))))
        for line in lines[2:end]]
end

epsilon_token(epsilon) = replace(@sprintf("%.4g", epsilon),
    "." => "p", "-" => "m")

function load_fold(epsilon)
    row = only(filter(r -> isapprox(r.epsilon,epsilon; atol=1e-13,rtol=0),
        numeric_csv(joinpath(SOURCE,"epsilon_fold_table.csv"))))
    op = BEKIsothermal._operators(DEGREE,2.0,0.6,0.5)
    profile = readdlm(joinpath(SOURCE,"profiles",
        "profile_epsilon_$(epsilon_token(epsilon)).csv"),',';skipstart=1)
    fields = permutedims(profile[:,2:5])
    solution = BEKConsistent.ConsistentSolution(row.Ro_f,
        2-row.Ro_f-row.Ro_f^2,1.0,0.72,op,fields,row.Tw_f,row.Hinf,
        row.residual,0)
    state = BEKConsistent.pack(fields)
    _,J = BEKConsistent._residual_jacobian(state,solution.Ro,solution.Tw,
        solution.gamma,solution.Pr,op)
    rownorm=max.([norm(view(J,i,:)) for i in axes(J,1)],1e-14)
    S=svd(J ./ rownorm)
    n=length(op.x)
    v=S.V[:,end]
    pivot=argmax(abs.(view(v,2n+1:3n)))+2n
    v[pivot]<0 && (v .*= -1)
    v ./= maximum(abs,view(v,2n+1:3n))
    w=S.U[:,end] ./ rownorm
    w ./= dot(w,v)
    seedfold=BEKConsistent.ConsistentFold(solution,v/norm(v),row.residual,0)
    corrected=BEKConsistent.solve_consistent_fold_at_h(-epsilon,seedfold;
        tolerance=2e-10)
    solution=corrected.solution
    state=BEKConsistent.pack(solution.fields)
    _,J=BEKConsistent._residual_jacobian(state,solution.Ro,solution.Tw,
        solution.gamma,solution.Pr,op)
    rownorm=max.([norm(view(J,i,:)) for i in axes(J,1)],1e-14)
    S=svd(J ./ rownorm)
    v=S.V[:,end]
    pivot=argmax(abs.(view(v,2n+1:3n)))+2n
    v[pivot]<0 && (v .*= -1)
    v ./= maximum(abs,view(v,2n+1:3n))
    w=S.U[:,end] ./ rownorm
    w ./= dot(w,v)
    solution,v,w
end

function branch_solution(fold,v,target_tw,seed_amplitude)
    fields=BEKConsistent.unpack(BEKConsistent.pack(fold.fields) .+
        seed_amplitude.*v,length(fold.operators.x))
    seed=BEKConsistent.ConsistentSolution(fold.Ro,fold.Co,fold.gamma,fold.Pr,
        fold.operators,fields,target_tw,fields[1,end],Inf,0)
    BEKConsistent.solve_consistent_fixed_tw(target_tw,seed;tolerance=5e-10)
end

function pencil_mode(solution,vfold,lambda0)
    J,M=BEKRadialResponse.radial_response_pencil(solution)
    q=copy(vfold)
    reference=vfold/dot(vfold,vfold)
    lambda=lambda0
    for iteration in 1:12
        A=J+lambda*M
        residual=vcat(A*q,dot(reference,q)-1)
        norm(residual,Inf)<2e-10 && return lambda,q,norm(A*q,Inf),iteration
        bordered=[A M*q; transpose(reference) 0.0]
        step=-(bordered\residual)
        q .+= view(step,1:length(q))
        lambda += step[end]
    end
    A=J+lambda*M
    lambda,q,norm(A*q,Inf),12
end

function main()
    rows=NamedTuple[]
    for (epsilon,target_A) in zip(EPSILONS,TARGET_AMPLITUDES)
        fold,v,w=load_fold(epsilon)
        state=BEKConsistent.pack(fold.fields)
        _,Mr=BEKRadialResponse.radial_response_pencil(fold)
        K=BEKConsistent._jacobian_directional_derivative(state,fold.Ro,
            fold.Tw,fold.gamma,fold.Pr,fold.operators,v)
        ctw=dot(w,BEKConsistent._tw_derivative(state,fold.operators))
        cr=dot(w,Mr*v)
        c2=0.5dot(w,K*v)
        delta_tw=-(c2/ctw)*target_A^2
        for signA in (-1.0,1.0)
            solution=branch_solution(fold,v,fold.Tw+delta_tw,signA*target_A)
            delta=BEKConsistent.pack(solution.fields)-state
            amplitude=dot(w,delta)
            predicted=-2c2*amplitude/cr
            lambda,q,pencil_residual,iterations=pencil_mode(solution,v,predicted)
            push!(rows,(N=DEGREE,epsilon,sign=Int(signA),Ro=fold.Ro,
                Tw_fold=fold.Tw,delta_Tw=delta_tw,amplitude,
                lambda_predicted=predicted,lambda_direct=lambda,
                relative_difference=abs(lambda-predicted)/max(abs(lambda),1e-30),
                branch_residual=solution.residual,pencil_residual,iterations))
            @printf("eps=%.3g sign=%+d A=%+.5e lambda=%+.5e pred=%+.5e rel=%.3e\n",
                epsilon,Int(signA),amplitude,lambda,predicted,
                rows[end].relative_difference)
        end
    end
    mkpath(dirname(OUTPUT))
    open(OUTPUT,"w") do io
        names=propertynames(first(rows)); println(io,join(string.(names),','))
        for row in rows
            for (i,name) in enumerate(names)
                i>1 && print(io,',')
                value=getproperty(row,name)
                value isa Integer ? print(io,value) : @printf(io,"%.12e",value)
            end
            println(io)
        end
    end
    println("Output: $OUTPUT")
end

main()

#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKRadialResponse.jl"))
using .BEKIsothermal, .BEKConsistent, .BEKRadialResponse
using LinearAlgebra
using Random
using Printf

function main()
    solution=BEKConsistent.solve_consistent_isothermal(-0.5;degree=60)
    D2,D1=BEKRadialResponse.radial_diffusion_operators(solution)
    n=length(solution.operators.x)
    rng=MersenneTwister(20260910)
    mode=randn(rng,4n)
    # Homogeneous boundary traces make replacement boundary rows vanish.
    for i in (1,n,n+1,2n,2n+1,3n,3n+1,4n)
        mode[i]=0
    end
    maximum_error=0.0
    for lambda in (-1.7,-0.3,0.0,0.8,2.0)
        analytic=(lambda^2*D2+lambda*D1)*mode
        direct=zeros(4n)
        rF=n+1:2n; rG=2n+1:3n; rT=3n+1:4n
        direct[rF] .= (lambda^2+2lambda) .* mode[rF]
        direct[rG] .= (lambda^2+2lambda) .* mode[rG]
        direct[rT] .= lambda^2 .* mode[rT]
        for row in (1,n,n+1,2n,2n+1,3n,3n+1,4n)
            direct[row]=0
        end
        error=norm(analytic-direct,Inf)
        maximum_error=max(maximum_error,error)
        @printf("lambda=%+.2f error=%.3e\n",lambda,error)
    end
    maximum_error <= 32eps(Float64) || error("radial diffusion operator failed")
end

main()

using Printf
include(joinpath(@__DIR__, "..", "src", "CompressibleBEKLinearStability.jl"))
using .CompressibleBEKLinearStability
refs=[("von_Karman",-1.0,0.38482,0.07759,285.36),("Ekman",0.0,0.2775,0.132,116.25),("Bodewadt",1.0,0.5231,0.127,27.26)]
for (name,Ro,a,b,R) in refs
    p=PhysicalParameters(Mr=0.3,Tw=1.0,Ro=Ro,R=R,beta=b,omega=0.0)
    n=NumericalParameters(N=59,domain=:rational,ymax=40.0,target=ComplexF64(a),neigs=12,tol=1e-10,maxit=1200,balance=true)
    c=prepare_case(p,n;allow_unverified_ro=(Ro != -1.0)); r=solve_mode(c)
    println(name)
    for j in eachindex(r.alpha)
        @printf("  %2d alpha=%+.9f%+.9fi |dref|=%.4e berr=%.2e\n",j,real(r.alpha[j]),imag(r.alpha[j]),abs(r.alpha[j]-a),r.backward_errors[j])
    end
end

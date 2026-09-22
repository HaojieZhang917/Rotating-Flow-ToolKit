using LinearAlgebra
using NonlinearEigenproblems
using Printf
include(raw"\\wsl.localhost\Ubuntu-24.04\home\zhj\Rotating-Flow-ToolKit\compress\BEK\roughness\CRD_STA.jl")

function old_case(Ro,Tw,Mr,R,beta,target)
    N=69; gamma=1.4; sigma=0.72; omega=0.0; Co=2-Ro-Ro^2
    u0,v0,w0,f,q,D,D2,x=baseflow_var(N,Ro,Co)
    H,T=T_ca(Mr,f,q,w0,gamma,Tw)
    F,G,H,T,rho,z=interp(u0,v0,H,T,x,N,"phy")
    lam=-(2/3).*T; kappa=T./sigma
    cof=Spatial_mode_BEK(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Mr/R,N,Ro,Co,D,D2)
    raw=assemble_mat(cof,D,D2,beta,omega)
    A=boudary_condition(raw...,N)
    vals,vecs=iar(PEP([A[1],A[2],A[3]]); σ=ComplexF64(target),neigs=6,maxit=800,tol=1e-10)
    return vals
end
for (name,Ro,a,b,R) in [("Ekman",0.0,0.2775,0.132,116.25),("Bodewadt",1.0,0.5231,0.127,27.26)]
    vals=old_case(Ro,1.0,0.3,R,b,a)
    println(name)
    for (j,v) in enumerate(vals)
        @printf(" %d %.12f%+.12fi |dref|=%.5e\n",j,real(v),imag(v),abs(v-a))
    end
end

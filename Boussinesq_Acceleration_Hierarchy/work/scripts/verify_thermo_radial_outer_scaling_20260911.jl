#!/usr/bin/env julia

using Random
using Printf

function full_residual(eps,X,Ro,Co,gamma,Pr,z)
    h,hY,f,fX,fY,fYY,fXX,g,gX,gY,gYY,gXX,t,tX,tY,tYY,tXX=z
    s=Co/2; omega=s+Ro; chi=1-gamma*(t-1); q=s+Ro*g
    H=eps*h; F=eps^2*f
    Heta=eps^2*hY
    Feta=eps^3*fY; Fetaeta=eps^4*fYY
    Fxi=eps^2*X*fX; Fxixi=eps^2*(X^2*fXX+X*fX)
    Geta=eps*gY; Getaeta=eps^2*gYY
    Gxi=X*gX; Gxixi=X^2*gXX+X*gX
    Teta=eps*tY; Tetaeta=eps^2*tYY
    Txi=X*tX; Txixi=X^2*tXX+X*tX
    kappa=eps^2/X^2
    (
      Heta+2*F+Fxi,
      Fetaeta+kappa*(Fxixi+2*Fxi)+
        chi*(Ro*(F^2+F*Fxi+H*Feta)-q^2/Ro)+omega^2/Ro,
      Getaeta+kappa*(Gxixi+2*Gxi)+
        chi*(Ro*(2*F*g+F*Gxi+H*Geta)+Co*F),
      Tetaeta+kappa*Txixi+Pr*Ro*(F*Txi+H*Teta),
    )
end

function scaled_outer_residual(eps,X,Ro,Co,gamma,Pr,z)
    h,hY,f,fX,fY,fYY,fXX,g,gX,gY,gYY,gXX,t,tX,tY,tYY,tXX=z
    s=Co/2; omega=s+Ro; chi=1-gamma*(t-1); q=s+Ro*g
    continuity=hY+2*f+X*fX
    algebraic=(-chi*q^2+omega^2)/Ro
    radial4=fYY+fXX+3*fX/X+chi*Ro*(f^2+X*f*fX+h*fY)
    azimuthal=gYY+gXX+3*gX/X+
      chi*(Ro*(2*f*g+X*f*gX+h*gY)+Co*f)
    thermal=tYY+tXX+tX/X+Pr*Ro*(X*f*tX+h*tY)
    (eps^2*continuity,algebraic+eps^4*radial4,
     eps^2*azimuthal,eps^2*thermal)
end

function main()
    rng=MersenneTwister(20260911)
    maximum_error=0.0
    for eps in (0.2,0.05,0.01), _ in 1:1000
        X=0.2+2rand(rng); Ro=-0.470367416464; Co=2-Ro-Ro^2
        z=randn(rng,17)
        a=full_residual(eps,X,Ro,Co,1.0,0.72,z)
        b=scaled_outer_residual(eps,X,Ro,Co,1.0,0.72,z)
        maximum_error=max(maximum_error,maximum(abs.(a.-b)))
    end
    @printf("thermo-radial outer scaling max identity error = %.3e\n",
        maximum_error)
    maximum_error<2e-12 || error("distinguished outer regression failed")
end

main()

#julia version:1.11.6
include("CRD_STA.jl")
using Plots
using LinearAlgebra
using NonlinearEigenproblems
using DelimitedFiles
using ProgressMeter
using PyCall
function ODE(N_cheb,Ro)
    Co = 2 - Ro - Ro^2
    u0,v0,w0,f,q,D,D2,x = baseflow_var(N_cheb,Ro,Co)
    return u0,v0,w0,f,q,D,D2,x
end
function BF(N_cheb,Tw,Mr,u0,v0,w0,f,q,D,D2,x)
    gamma = 1.4
    sigma = 0.72
    H,T = T_ca(Mr,f,q,w0,gamma,Tw)
    F,G,H,T,rho,z = interp(u0,v0,H,T,x,N_cheb,"phy")
    lam = - (2/3) * T
    kappa = (1/sigma) * T
    return F,G,H,T,rho,z,lam,kappa,D,D2
end
N_cheb = 199
Ro = -1.0
Tw = 1.04
Mr = 0.3
gamma = 1.4
sigma = 0.72
Co = 2-Ro-Ro^2
be = 0.07759
num = 1
omega = 0.0
data = [0 0 0 0 0]
R = 285.36
Ma = Mr/R
u0,v0,w0,f,q,D,D2,x = ODE(N_cheb,Ro)
F,G,H,T,rho,z,lam,kappa,D,D2 = BF(N_cheb,Tw,Mr,u0,v0,w0,f,q,D,D2,x)
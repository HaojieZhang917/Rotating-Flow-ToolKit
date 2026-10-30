include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra

const OUTDIR=joinpath(@__DIR__,"results","wavenumber_matching")
const R0=285.36; const BETA=0.07759; const OMEGA=-BETA; const TARGET=0.18482
const AINC=1e-3+0im; const EPSROUGH=1e-3; const ZREF=10.0; const ZI=10.0; const N=109; const L=30.0
make_p(a)=LSTParameters(mode=:spatial,R=R0,alpha=a,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im)
function incident(base,p)
    ms=farfield_six_state_modes(base;p=p); cand=[m for m in ms if imag(m.kz)<-1e-8]; low=minimum(abs(imag(x.kz)) for x in cand); near=[x for x in cand if abs(abs(imag(x.kz))-low)<1e-8]; m=near[argmax([real(x.kz) for x in near])]
    y=copy(m.state); y./=norm(y[1:3]); AINC.*y.*exp(im*m.kz*(L-ZREF))
end
function profile4(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z); out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n]); for j in eachindex(zq); out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j]); end
    end; out
end
function energy(resp); zq,w=subinterval_cc_grid(ZI,128); q=profile4(resp,zq); sum(w.*vec(sum(abs2.(q[1:3,:]);dims=1))); end
base=solve_von_karman_baseflow(); alphas=[.25,.30,.35,.38482,.45]; io=IOBuffer(); println(io,"alpha_f,kappa_abs,kappa_signed,target_alpha,Gv,Gw,Gf,H_eff,source_norm,rho")
for a in alphas
    k_abs=a-TARGET; k_signed=-abs(k_abs); @assert isapprox(a+k_signed,TARGET;atol=1e-12)
    p=make_p(a); yi=incident(base,p); s=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yi); qf=as_four_state_response(s)
    qh=solve_roughness_distortion(base;p=p,kappa=k_signed,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
    co=solve_coupled_response(base,qh.response,qf;roughness_amplitude=EPSROUGH,zmax=L)
    Ev=energy(co.volume_response); Ew=energy(co.wall_response); Ef=energy(co.response); Ec=Ef-Ev-Ew; src=norm(co.source[1:3*length(co.response.z)])
    he=sqrt(Ev)/src; rho=Ew>0 ? Ec/(2sqrt(Ev*Ew)) : NaN
    println(io,join([a,k_abs,k_signed,TARGET,Ev/(EPSROUGH^2*abs(AINC)^2),Ew/(EPSROUGH^2*abs(AINC)^2),Ef/(EPSROUGH^2*abs(AINC)^2),he,src,rho],","))
end
write(joinpath(OUTDIR,"fixed_target_matching.csv"),String(take!(io))); println("wrote fixed target")

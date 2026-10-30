include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra

const OUTDIR=joinpath(@__DIR__,"results","roughness_wavenumber_scan")
const R0=285.36; const ALPHA=0.38482; const BETA=0.07759; const OMEGA=-BETA
const ZREF=10.0; const AINC=1e-3+0im; const EPSROUGH=1e-3; const N=109; const L=30.0
const KLIST=[0.05,0.4,1.5]

function incident_state(base,p,z)
    modes=farfield_six_state_modes(base;p=p); m=modes[argmin([abs(x.kz+0.1491006263im) for x in modes])]
    y=copy(m.state); y./=norm(y[1:3]); return AINC.*y.*exp(im*m.kz*(z-ZREF))
end

function profile(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z)
    out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n])
        for j in eachindex(zq); out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j]); end
    end
    out
end

base=solve_von_karman_baseflow()
p=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=0.38+0im)
yinc=incident_state(base,p,L)
d6=solve_six_state_lns(base;p=p,zmax=L,boundary=:dirichlet,incident_state=yinc,edge_velocity=yinc[1:3])
m6=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yinc)
qfd=as_four_state_response(d6); qfm=as_four_state_response(m6)
zq=collect(range(0.0,10.0;length=160))
io=IOBuffer(); println(io,"kappa,z,uD_re,uD_im,vD_re,vD_im,wD_re,wD_im,uM_re,uM_im,vM_re,vM_im,wM_re,wM_im")
for k in KLIST
    qh=solve_roughness_distortion(base;p=p,kappa=k,n=0.0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
    d=solve_coupled_response(base,qh.response,qfd;roughness_amplitude=EPSROUGH,zmax=L)
    m=solve_coupled_response(base,qh.response,qfm;roughness_amplitude=EPSROUGH,zmax=L)
    pd=profile(d.response,zq); pm=profile(m.response,zq)
    for j in eachindex(zq)
        println(io,join([k,zq[j],real(pd[1,j]),imag(pd[1,j]),real(pd[2,j]),imag(pd[2,j]),real(pd[3,j]),imag(pd[3,j]),
                         real(pm[1,j]),imag(pm[1,j]),real(pm[2,j]),imag(pm[2,j]),real(pm[3,j]),imag(pm[3,j])],","))
    end
end
write(joinpath(OUTDIR,"roughness_scan_profiles.csv"),String(take!(io)))
println("wrote profiles")

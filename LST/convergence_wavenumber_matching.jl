include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra

const OUTDIR=joinpath(@__DIR__,"results","wavenumber_matching")
const R0=285.36; const BETA=0.07759; const OMEGA=-BETA; const AINC=1e-3+0im; const EPSROUGH=1e-3; const ZREF=10.0; const ZI=10.0
make_p(a,N)=LSTParameters(mode=:spatial,R=R0,alpha=a,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im)
function incident(base,p,L)
    ms=farfield_six_state_modes(base;p=p); cand=[m for m in ms if imag(m.kz)<-1e-8]
    lowabs=minimum(abs(imag(x.kz)) for x in cand); near=[x for x in cand if abs(abs(imag(x.kz))-lowabs)<1e-8]
    m=near[argmax([real(x.kz) for x in near])]
    y=copy(m.state); y./=norm(y[1:3]); AINC.*y.*exp(im*m.kz*(L-ZREF))
end
function profile4(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z); out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n]); for j in eachindex(zq); out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j]); end
    end; out
end
function energy(resp,L)
    zq,w=subinterval_cc_grid(L,128); q=profile4(resp,zq); sum(w.*vec(sum(abs2.(q[1:3,:]);dims=1)))
end
function run(base,a,k,N,L,branch)
    p=make_p(a,N); yi=incident(base,p,L); s=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yi); qf=as_four_state_response(s)
    ks=branch==:plus ? k : -k; qh=solve_roughness_distortion(base;p=p,kappa=ks,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
    co=solve_coupled_response(base,qh.response,qf;roughness_amplitude=EPSROUGH,zmax=L)
    return (alpha_f=a,kappa=k,branch=String(branch),N=N,zmax=L,Gf=energy(co.response,ZI)/(EPSROUGH^2*abs(AINC)^2),Gv=energy(co.volume_response,ZI)/(EPSROUGH^2*abs(AINC)^2),Gw=energy(co.wall_response,ZI)/(EPSROUGH^2*abs(AINC)^2),qf_deriv=qf.velocity_derivative_difference,residual=norm(co.response.interior_residual,Inf))
end
base=solve_von_karman_baseflow(); tests=[(.38482,.005,:plus),(.38482,.4,:plus),(.38482,.2,:minus),(.2,.02,:plus),(.45,.02,:plus)]
rows=NamedTuple[]
for (a,k,b) in tests, N in [89,109,129], L in [20.0,30.0,40.0]
    println("a=",a," k=",k," b=",b," N=",N," L=",L); push!(rows,run(base,a,k,N,L,b))
end
cols=String.(fieldnames(typeof(rows[1]))); io=IOBuffer(); println(io,join(cols,",")); for r in rows; println(io,join([string(getproperty(r,Symbol(c))) for c in cols],",")); end
write(joinpath(OUTDIR,"wavenumber_matching_convergence.csv"),String(take!(io))); println("wrote convergence")

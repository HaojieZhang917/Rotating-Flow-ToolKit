include(joinpath(@__DIR__,"VonKarmanLST.jl")); using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl")); using .SixStateLNS
using LinearAlgebra
const OUTDIR=joinpath(@__DIR__,"results","second_order_modal"); mkpath(OUTDIR)
const R0=285.36; const ALPHA=.38482; const BETA=.07759; const OMEGA=-BETA; const AINC=1e-3+0im; const EPS=1e-3; const ZI=10.0; const K=.185
function inc(base,p,L)
 ms=farfield_six_state_modes(base;p=p); c=[m for m in ms if imag(m.kz)<-1e-8]; lm=minimum(abs(imag(m.kz)) for m in c); c=[m for m in c if abs(abs(imag(m.kz))-lm)<1e-8]; m=c[argmax([real(x.kz) for x in c])]; y=m.state/norm(m.state[1:3]); AINC*y*exp(im*m.kz*(L-10.0))
end
function e(resp)
 z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z); x,w=subinterval_cc_grid(ZI,128); q=zeros(ComplexF64,3,length(x));
 for c in 1:3
  val=resp.q[(c-1)*n+1:c*n]; for j in eachindex(x); q[c,j]=VonKarmanLST._barycentric_interpolate(z,val,bw,x[j]); end
 end
 sum(w.*vec(sum(abs2.(q);dims=1)))
end
base=solve_von_karman_baseflow(); rows=NamedTuple[]
for N in [109,129,169], L in [20.0,30.0,40.0,50.0]
 p=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im); println("N=",N," L=",L)
 qf=as_four_state_response(solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=inc(base,p,L)))
 qh=solve_roughness_distortion(base;p=p,kappa=-K,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
 s=solve_six_state_coupled_response(base,qh.response,qf;roughness_amplitude=EPS,zmax=L,boundary=:modal)
 push!(rows,(N=N,zmax=L,G= e(as_four_state_response(s.response))/(EPS^2*abs(AINC)^2),bc=norm(s.response.modal.C*s.response.outer_state),prim=s.response.residuals.primitive_inf))
end
cols=String.(fieldnames(typeof(rows[1]))); io=IOBuffer(); println(io,join(cols,",")); for r in rows; println(io,join([string(getproperty(r,Symbol(c))) for c in cols],",")); end; write(joinpath(OUTDIR,"second_order_modal_convergence.csv"),String(take!(io)))

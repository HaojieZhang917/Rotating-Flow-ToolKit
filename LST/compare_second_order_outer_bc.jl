include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra

const OUTDIR=joinpath(@__DIR__,"results","second_order_modal")
mkpath(OUTDIR)
const R0=285.36; const ALPHA=0.38482; const BETA=0.07759; const OMEGA=-BETA
const AINC=1e-3+0im; const EPSROUGH=1e-3; const ZREF=10.0; const ZI=10.0; const N=109; const L=30.0
make_p()=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im)
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
function energy(resp)
    zq,w=subinterval_cc_grid(ZI,128); q=profile4(resp,zq); sum(w.*vec(sum(abs2.(q[1:3,:]);dims=1)))
end
function metrics(resp)
    E=energy(resp); (G=E/(EPSROUGH^2*abs(AINC)^2),E=E)
end
function dynamic_interior(sol)
    n=length(sol.z); rows=setdiff(1:length(sol.Y),[1,n+1,2n+1,n,2n,3n]); norm(sol.residuals.dynamic[rows],Inf)
end
function run_k(base,cache,k)
    qh=solve_roughness_distortion(base;p=cache.p,kappa=-k,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
    old=solve_coupled_response(base,qh.response,cache.qf;roughness_amplitude=EPSROUGH,zmax=L)
    sixD=solve_six_state_coupled_response(base,qh.response,cache.qf;roughness_amplitude=EPSROUGH,zmax=L,boundary=:dirichlet)
    sixM=solve_six_state_coupled_response(base,qh.response,cache.qf;roughness_amplitude=EPSROUGH,zmax=L,boundary=:modal)
    od=metrics(old.response); d=metrics(as_four_state_response(sixD.response)); m=metrics(as_four_state_response(sixM.response))
    zq=collect(range(0.0,ZI;length=128)); po=profile4(old.response,zq); pd=profile4(as_four_state_response(sixD.response),zq); pm=profile4(as_four_state_response(sixM.response),zq)
    return (kappa=k,G_old=od.G,G_D=d.G,G_M=m.G,old_D_profile=norm(po-pd)/max(norm(po),eps(Float64)),D_M_profile=norm(pd-pm)/max(norm(pm),eps(Float64)),
            dynamic_D=dynamic_interior(sixD.response),dynamic_M=dynamic_interior(sixM.response),primitive_D=sixD.response.residuals.primitive_inf,primitive_M=sixM.response.residuals.primitive_inf,
            modal_bc=norm(sixM.response.modal.C*sixM.response.outer_state),qh_residual=norm(qh.response.interior_residual,Inf),qh_bc=norm(qh.response.boundary_error,Inf))
end
base=solve_von_karman_baseflow(); p=make_p(); yi=incident(base,p); s=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yi); cache=(p=p,qf=as_four_state_response(s))
rows=NamedTuple[]
for k in [0.15,0.18,0.185,0.19,0.20,0.25]
    println("kappa=",k); push!(rows,run_k(base,cache,k))
end
cols=String.(fieldnames(typeof(rows[1]))); io=IOBuffer(); println(io,join(cols,",")); for r in rows; println(io,join([string(getproperty(r,Symbol(c))) for c in cols],",")); end
write(joinpath(OUTDIR,"second_order_modal_comparison.csv"),String(take!(io)))
println("wrote second-order comparison")

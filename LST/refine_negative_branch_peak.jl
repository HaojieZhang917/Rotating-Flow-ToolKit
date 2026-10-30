include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra

const OUTDIR=joinpath(@__DIR__,"results","negative_branch_peak")
mkpath(OUTDIR)
const R0=285.36; const ALPHA=0.38482; const BETA=0.07759; const OMEGA=-BETA
const AINC=1e-3+0im; const EPSROUGH=1e-3; const ZREF=10.0; const ZI=10.0
const N=109; const L=30.0; const KLIST=collect(0.14:0.005:0.28)

make_p()=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im)
function incident(base,p)
    ms=farfield_six_state_modes(base;p=p); cand=[m for m in ms if imag(m.kz)<-1e-8]
    low=minimum(abs(imag(x.kz)) for x in cand); near=[x for x in cand if abs(abs(imag(x.kz))-low)<1e-8]
    m=near[argmax([real(x.kz) for x in near])]; y=copy(m.state); y./=norm(y[1:3])
    AINC.*y.*exp(im*m.kz*(L-ZREF))
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
function inner4(a,b)
    zq,w=subinterval_cc_grid(ZI,128); pa=profile4(a,zq); pb=profile4(b,zq); sum(w.*vec(sum(conj.(pa[1:3,:]).*pb[1:3,:];dims=1)))
end
function run(base,cache,k)
    ks=-k; qh=solve_roughness_distortion(base;p=cache.p,kappa=ks,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
    co=solve_coupled_response(base,qh.response,cache.qf;roughness_amplitude=EPSROUGH,zmax=L)
    vol=co.volume_response; wall=co.wall_response; full=co.response
    Ev=energy(vol); Ew=energy(wall); Ef=energy(full); Ec=2real(inner4(vol,wall)); src=norm(co.source[1:3*length(co.response.z)])
    return (kappa=k,target_alpha=ALPHA-k,Gv=Ev/(EPSROUGH^2*abs(AINC)^2),Gw=Ew/(EPSROUGH^2*abs(AINC)^2),Gf=Ef/(EPSROUGH^2*abs(AINC)^2),
            Ev=Ev,Ew=Ew,Ef=Ef,Ec=Ec,rho=Ew>0 ? Ec/(2sqrt(Ev*Ew)) : NaN,Rwall=Ew/Ev,
            source_norm=src,H_eff=sqrt(Ev)/src,linear_error=norm(profile4(full,collect(range(0.0,ZI;length=100)))-profile4(vol,collect(range(0.0,ZI;length=100)))-profile4(wall,collect(range(0.0,ZI;length=100))))/sqrt(Ef),
            residual=norm(co.response.interior_residual,Inf),qh_residual=norm(qh.response.interior_residual,Inf),qh_bc=norm(qh.response.boundary_error,Inf),
            co=co,qh=qh)
end
base=solve_von_karman_baseflow(); p=make_p(); yi=incident(base,p); sol=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yi); cache=(p=p,qf=as_four_state_response(sol))
rows=NamedTuple[]
for k in KLIST
    println("kappa=",k); r=run(base,cache,k)
    push!(rows,(kappa=r.kappa,target_alpha=r.target_alpha,Gv=r.Gv,Gw=r.Gw,Gf=r.Gf,Ev=r.Ev,Ew=r.Ew,Ef=r.Ef,Ec=r.Ec,rho=r.rho,Rwall=r.Rwall,source_norm=r.source_norm,H_eff=r.H_eff,linear_error=r.linear_error,residual=r.residual,qh_residual=r.qh_residual,qh_bc=r.qh_bc))
end
cols=String.(fieldnames(typeof(rows[1]))); io=IOBuffer(); println(io,join(cols,",")); for r in rows; println(io,join([string(getproperty(r,Symbol(c))) for c in cols],",")); end
write(joinpath(OUTDIR,"negative_branch_refined.csv"),String(take!(io)))
println("wrote refined scan")

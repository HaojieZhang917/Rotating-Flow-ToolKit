include(joinpath(@__DIR__,"VonKarmanLST.jl")); using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl")); using .SixStateLNS
using Serialization, LinearAlgebra
x=deserialize(joinpath(@__DIR__,"results","roughness_reconstruction","reconstructed_response.bin"))
base=solve_von_karman_baseflow(); p=LSTParameters(mode=:spatial,R=285.36,alpha=.38482,beta=.07759,omega=-.07759,N=x.params.N,nev=8,spatial_shift=.38+0im)
ms=farfield_six_state_modes(base;p=p); c=[m for m in ms if imag(m.kz)<-1e-8]; lm=minimum(abs(imag(m.kz)) for m in c); c=[m for m in c if abs(abs(imag(m.kz))-lm)<1e-8]; m=c[argmax([real(q.kz) for q in c])]
yi=1e-3*m.state/norm(m.state[1:3])*exp(im*m.kz*(x.params.zmax-10)); qf=as_four_state_response(solve_six_state_lns(base;p=p,zmax=x.params.zmax,boundary=:modal,incident_state=yi))
err=Float64[]
for i in eachindex(x.r)
    h=gaussian_roughness(x.r[i];rf=0,c=1,epsilon=1e-3,branch=:real); ref=-h.*qf.wall_shear.*exp(im*.38482*x.r[i]); push!(err,norm(x.q[1:3,i,1]-ref)/max(norm(ref),eps(Float64)))
end
println(maximum(err))

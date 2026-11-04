include(joinpath(@__DIR__,"VonKarmanLST.jl")); using .VonKarmanLST
base=solve_von_karman_baseflow(); rows=[]
for beta in [.04,.06,.07759,.10,.12], omega in [-.05,0.0,.05]
 p=LSTParameters(mode=:spatial,R=285.36,alpha=.38482,beta=beta,omega=omega,N=39,nev=4*40,spatial_shift=.38+0im); b=assemble_spatial_blocks(base;p=p,grid=:finite,zmax=30.0); L0,L1,L2=VonKarmanLST._spatial_pencil(b,p); k=b.keep; op=(L0=Matrix{ComplexF64}(L0[k,k]),L1=Matrix{ComplexF64}(L1[k,k]),L2=Matrix{ComplexF64}(L2[k,k]),z=b.z,D=b.D,params=p,blocks=b); m=solve_spatial(op;shift=.38+0im,nev=size(op.L0,1)); c=[x for x in m if real(x.alpha)>.2&&real(x.alpha)<.7]; isempty(c) && continue; x=c[argmin([abs(q.alpha-.38482) for q in c])]; push!(rows,(beta=beta,omega=omega,ar=real(x.alpha),ai=imag(x.alpha)))
end
io=IOBuffer(); println(io,"beta,omega,alpha_real,alpha_imag"); for r in rows; println(io,"$(r.beta),$(r.omega),$(r.ar),$(r.ai)"); end; write(joinpath(@__DIR__,"results","spatial_lstreceptivity","candidate_parameter_scan.csv"),String(take!(io)))

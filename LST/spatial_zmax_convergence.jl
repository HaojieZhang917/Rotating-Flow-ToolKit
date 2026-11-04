include(joinpath(@__DIR__,"VonKarmanLST.jl")); using .VonKarmanLST; using LinearAlgebra
const R=285.36; const beta=.07759; const omega=0.0; const N=69
base=solve_von_karman_baseflow(); rows=[]
for zmax in [20.0,30.0,40.0,50.0]
 p=LSTParameters(mode=:spatial,R=R,alpha=.38482,beta=beta,omega=omega,N=N,nev=4*(N+1),spatial_shift=.38+0im)
 b=assemble_spatial_blocks(base;p=p,grid=:finite,zmax=zmax); L0,L1,L2=VonKarmanLST._spatial_pencil(b,p); k=b.keep; op=(L0=Matrix{ComplexF64}(L0[k,k]),L1=Matrix{ComplexF64}(L1[k,k]),L2=Matrix{ComplexF64}(L2[k,k]),z=b.z,D=b.D,params=p,blocks=b,zmax=zmax)
 modes=solve_spatial(op;shift=.38+0im,nev=size(op.L0,1)); c=[m for m in modes if real(m.alpha)>.2&&real(m.alpha)<.7]; m=c[argmin([abs(x.alpha-.38482) for x in c])]; res=norm((op.L0+m.alpha*op.L1+m.alpha^2*op.L2)*m.vector)/max(norm(op.L0)*norm(m.vector),eps(Float64)); push!(rows,(zmax=zmax,alpha_real=real(m.alpha),alpha_imag=imag(m.alpha),residual=res))
end
io=IOBuffer(); println(io,"zmax,alpha_real,alpha_imag,residual"); for r in rows; println(io,"$(r.zmax),$(r.alpha_real),$(r.alpha_imag),$(r.residual)"); end; write(joinpath(@__DIR__,"results","spatial_lstreceptivity","spatial_zmax_convergence.csv"),String(take!(io)))

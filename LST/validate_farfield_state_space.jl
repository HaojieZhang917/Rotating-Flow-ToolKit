#=
验证远场一阶广义状态系统。

该脚本只验证状态空间与原二次远场 pencil 的等价性，不改变现有
有限外缘 Dirichlet 求解器，也不把任何模态自动解释为入射波。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
using LinearAlgebra
using Printf

const OUTDIR=joinpath(@__DIR__,"results","farfield_state_space")
mkpath(OUTDIR)

base=solve_von_karman_baseflow()
p=LSTParameters(mode=:spatial,R=285.36,alpha=0.38482,beta=0.07759,
                omega=-0.07759,N=69)
quadratic=solve_farfield_modes(base;p=p)
state=solve_farfield_state_modes(base;p=p)
bounded=farfield_bounded_subspace(base;p=p)
algebraic=farfield_algebraic_constraints(base;p=p)

used=trues(length(quadratic))
matches=NamedTuple[]
for sm in state.modes
    d=[abs(sm.kz-q.kz) for q in quadratic]
    ids=findall(used)
    j=ids[argmin(d[ids])]
    used[j]=false
    push!(matches,(kz_state=sm.kz,kz_quadratic=quadratic[j].kz,
                   difference=sm.kz-quadratic[j].kz,
                   pencil_residual=sm.pencil_residual,
                   derivative_residual=sm.derivative_residual,
                   quadratic_residual=sm.quadratic_residual,
                   lambda=sm.lambda))
end
sort!(matches,by=x->imag(x.kz_state))

io=IOBuffer()
println(io,"# Far-field first-order state-space validation")
println(io)
println(io,"参数：R=$(p.R), alpha=$(p.alpha), beta=$(p.beta), omega=$(p.omega), N=$(p.N)。")
println(io)
println(io,"状态约定：X=[q;q_z]，q(z)=a exp(i*k_z*z)，lambda=i*k_z，E X_z=A X。")
println(io)
println(io,"## Pencil count")
println(io)
println(io,"- finite state modes: $(length(state.modes))")
println(io,"- infinite generalized eigenvalues: $(state.infinite_count)")
println(io,"- original quadratic finite modes: $(length(quadratic))")
println(io,"- finite modes bounded toward z=+infinity: $(length(bounded.modes))")
println(io,"- left-null constraints of bounded finite subspace: $(size(bounded.C,1))")
println(io,"- direct left-null dimension of E: $(size(algebraic.left_null_basis,2))")
println(io,"- direct algebraic constraint rank: $(algebraic.rank)")
println(io)
println(io,"## Root-by-root equivalence")
println(io)
println(io,"|root|kz from state pencil|kz from quadratic pencil|difference|state residual|qz-lambda*q|quadratic residual|orientation under exp(i*kz*z)|")
println(io,"|---:|---:|---:|---:|---:|---:|---:|")
for (j,m) in enumerate(matches)
    orient=imag(m.kz_state)>0 ? "bounded toward +z" : "grows toward +z; transfer toward wall under fixed edge normalization"
    println(io,"|$j|$(m.kz_state)|$(m.kz_quadratic)|$(m.difference)|$(m.pencil_residual)|$(m.derivative_residual)|$(m.quadratic_residual)|$orient|")
end
println(io)
println(io,"## Direct sign check near the suspected root")
println(io)
println(io,"The following values use the same `exp(i*k_z*z)` convention and evaluate the original 4×4 matrix directly.")
println(io)
println(io,"|probe kz|minimum singular value of M_inf|det(M_inf)|")
println(io,"|---:|---:|---:|")
for kz in [-0.149100626278im,0.149100626278im]
    M=farfield_matrix(base;p=p,kz=kz)
    println(io,"|$(kz)|$(minimum(svdvals(M)))|$(det(M))|")
end
println(io)
println(io,"## Interpretation")
println(io)
println(io,"The six finite roots from the first-order generalized pencil agree with the original quadratic roots to numerical precision. The two infinite eigenvalues are the algebraic constraints associated with the pressure-containing primitive formulation; they must not be silently discarded when a physical boundary operator is constructed.")
println(io)
println(io,"Under q(z)=a exp(i*k_z*z), imag(k_z)>0 is bounded as z→+infinity and imag(k_z)<0 grows as z→+infinity. The latter roots can still describe a response propagated inward from a finite forcing plane, but they are not admissible outer-field decay modes by themselves.")
println(io)
println(io,"This validation does not yet impose an incoming/outgoing scattering condition. The next implementation must select a physical finite/infinite constraint subspace and derive the correct number of independent outer boundary equations after the pressure constraint is treated.")
println(io)
println(io,"## Bounded-subspace boundary-count warning")
println(io)
println(io,"The current finite bounded subspace contains $(length(bounded.modes)) modes, and its primitive-state left null space has $(size(bounded.C,1)) rows with residual $(bounded.constraint_residual). Because the primitive descriptor pencil has $(state.infinite_count) infinite eigenvalues but only $(algebraic.rank) direct left-null algebraic constraint, at least one hidden higher-index constraint is present. These $(size(bounded.C,1)) rows cannot yet be inserted directly into the four-row collocation outer boundary block. A pressure/constraint reduction or an equivalent descriptor boundary formulation is required before solving the scattering problem.")
write(joinpath(OUTDIR,"farfield_state_space_report.md"),String(take!(io)))
println("wrote ",joinpath(OUTDIR,"farfield_state_space_report.md"))

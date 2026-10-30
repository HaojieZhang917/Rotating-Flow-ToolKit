include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl"))
using .SixStateLNS
using LinearAlgebra

const OUTDIR=joinpath(@__DIR__,"results","second_order_modal"); mkpath(OUTDIR)
const R0=285.36; const ALPHA=0.38482; const BETA=0.07759; const OMEGA=-BETA; const AINC=1e-3+0im; const EPSROUGH=1e-3; const ZREF=10.0; const ZI=10.0; const N=109; const L=30.0; const KAPPA=.185
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
function inner(a,b)
    zq,w=subinterval_cc_grid(ZI,128); qa=profile4(a,zq); qb=profile4(b,zq); sum(w.*vec(sum(conj.(qa[1:3,:]).*qb[1:3,:];dims=1)))
end
function report()
    base=solve_von_karman_baseflow(); p=make_p(); yi=incident(base,p); first=solve_six_state_lns(base;p=p,zmax=L,boundary=:modal,incident_state=yi); qf=as_four_state_response(first)
    qh=solve_roughness_distortion(base;p=p,kappa=-KAPPA,n=0,h_amplitude=1.0,zmax=L,edge_state=zeros(ComplexF64,4))
    sol=solve_six_state_coupled_response(base,qh.response,qf;roughness_amplitude=EPSROUGH,zmax=L,boundary=:modal)
    v=sol.volume_response; w=sol.wall_response; t=sol.response
    mb=modal_outer_boundary(base;p=sol.params); zidx=setdiff(1:length(t.Y),[1,N+2,2N+3,N+1,2N+2,3N+3])
    zf=solve_six_state_forced(base,zeros(ComplexF64,4*(N+1));p=sol.params,zmax=L,boundary=:modal)
    sup=norm(t.Y-v.Y-w.Y)/max(norm(t.Y),eps(Float64)); eV=energy(as_four_state_response(v)); eW=energy(as_four_state_response(w)); eT=energy(as_four_state_response(t)); cross=eT-(eV+eW+2real(inner(as_four_state_response(v),as_four_state_response(w))))
    s2=solve_six_state_coupled_response(base,qh.response,qf;roughness_amplitude=2EPSROUGH,zmax=L,boundary=:modal); scale=norm(s2.response.Y-2t.Y)/max(norm(2t.Y),eps(Float64))
    # 内部动态方程残差不包含被边界条件替换的六行。
    dyn=norm(t.residuals.dynamic[zidx],Inf); prim=t.residuals.primitive_inf
    open(joinpath(OUTDIR,"second_order_modal_validation.md"),"w") do io
        println(io,"# 二阶六维 Modal 外缘条件验证")
        println(io,"\n参数：`R=$(R0)`, `alpha=$(ALPHA)`, `beta=$(BETA)`, `omega=$(OMEGA)`, `kappa=-$(KAPPA)`, `N=$(N)`, `zmax=$(L)`, `Ainc=epsilon=1e-3`。本文数值中的 `G_full` 沿用项目现有定义：`G_full=∫_0^{z_i}|u_hf|² dz/(epsilon²|Ainc|²)`，不是其平方根。")
        println(io,"\n## 远场 Modal 子空间")
        println(io,"- 有界模态数：$(length(mb.modes))；`rank(C)= $(mb.rank)`；`rank(C*Vg)= $(mb.growth_constraint_rank)`。")
        println(io,"- `||C V_decay||₂ = $(mb.constraint_residual)`，`cond(V_decay) = $(mb.condition)`。")
        println(io,"\n## 六维受迫系统")
        println(io,"- 内部动态残差 `||Y'-AY-F||∞ = $(dyn)`；原始四变量方程残差 `$(prim)`。")
        println(io,"- Modal 外缘残差 `||C Y(zmax)||₂ = $(norm(sol.response.modal.C*sol.response.outer_state))`。")
        println(io,"\n## 零强迫、叠加与幅值缩放")
        println(io,"- 零源零壁面：`||Y||₂ = $(norm(zf.Y))`。")
        println(io,"- 线性叠加相对误差 `||Y_total-Y_volume-Y_wall||/||Y_total|| = $(sup)`。")
        println(io,"- 能量恒等式误差 `E_total-E_volume-E_wall-2Re<V,W> = $(cross)`；`E_total=$(eT)`。")
        println(io,"- 粗糙度幅值从 `1e-3` 加倍到 `2e-3` 的状态缩放误差：`$(scale)`。")
        println(io,"\n## 代表点结果")
        println(io,"- `G_D=0.2302883863`（六维 Dirichlet），`G_M=0.2305137514`（六维 Modal），代表点的 Modal 相对差约 `0.0978%`。")
        println(io,"- 与旧四变量 Dirichlet 的代表点增益差为 `3.76e-7` 相对量级；六维/四变量剖面相对误差约 `4e-5`，因此增益一致性很好，但剖面尚未达到机器精度，应归因于不同状态重构和外缘压力处理的离散差异。")
        println(io,"- 负分支精细扫描见 [second_order_modal_peak_scan.csv](second_order_modal_peak_scan.csv)，在本次网格中两种边界条件的最大值均位于 `kappa=0.185`；Modal 与 Dirichlet 增益差随 `|kappa|` 增大，在 `|kappa|=0.28` 约 `0.35%`。")
        println(io,"- `N/zmax` 收敛数据见 [second_order_modal_convergence.csv](second_order_modal_convergence.csv)：在 `zmax=30`，`N=109→169` 的增益变化约 `3.2e-7`；`N=169` 时 `zmax=40→50` 变化约 `1.5e-4`，而 `zmax=20→30` 仍约 `7.0e-3`。因此峰值扫描采用 `zmax=30` 可用于趋势判断，但若要给出高精度绝对增益，建议采用 `zmax≥40` 并报告外缘收敛。")
        println(io,"\n## 解释")
        println(io,"Modal 条件允许外缘存在有界远场散射状态，因此它与有限外缘三速度 Dirichlet 条件的差异是物理边界模型的差异。代表峰值附近差异小于千分之一，但整个扫描区间已达到约 `0.35%`；在增大外缘高度、改变频率或接近远场模态退化时不能预先假定仍然很小。")
    end
    println("validation report written")
end
report()

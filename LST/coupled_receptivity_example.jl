#=
有限外缘 Dirichlet + Gaussian 粗糙度的 frozen-R 二阶耦合示例。

该示例依次计算 q_f、单位粗糙度 q_h、交叉体积源 F_hf 和 q_hf，
并用纯 Julia 写出一个 SVG 结果图。它不代表真实自由涡波入射或
全局径向非平行旋转盘解；Re_h、r_f、c 和周向 beta 映射必须按具体
论文参数替换。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
using LinearAlgebra

# ---------------------------- 参数 ----------------------------
R=285.36
N=29                         # 先用较低阶数做串联验证；收敛计算改为 69 或更高
zmax=20.0
alpha_f=0.38482
beta_f=0.07759
omega_f=-0.07759
edge_amplitude=1e-3

# 这些是示例值，不能替代论文中的 Re_h、r_f、c 定义。
Re_h=1000.0
epsilon=1.0/Re_h
rf=0.0
c=0.5
n_rough=0                         # n=5 需要先提供已核对的 beta 映射
beta_h=0.0                        # n=0 时没有周向波数换算问题
nk=32                         # 取偶数以避免 κ=0 的未规范定常零模

base=solve_von_karman_baseflow()
p_f=LSTParameters(mode=:spatial,R=R,alpha=alpha_f,beta=beta_f,
                  omega=omega_f,N=N)
kvec=[alpha_f,beta_f,0.0]
edge_state=first(vortical_edge_states(kvec;amplitude=edge_amplitude))
qf=solve_forced_spatial(base,edge_state;p=p_f,zmax=zmax)

# Gaussian 频谱使用 Hhat=∫Hr exp(-im*kappa*r)dr，逆变换带 1/(2π)。
grid=roughness_fourier_grid(c=c,kappa_max=6.0/c,nk=nk)
qh_modes=Any[]
qhf_modes=Any[]
qv_modes=Any[]
qw_modes=Any[]
source_modes=Any[]

for kappa in grid.kappa
    rough=solve_roughness_distortion(base;p=p_f,kappa=kappa,n=n_rough,
                                     beta_shift=beta_h,h_amplitude=1.0,
                                     zmax=zmax)
    push!(qh_modes,rough.response)
    # 这里求单位 Fourier 粗糙度系数；epsilon*Hhat*Δkappa/(2π)
    # 只在最后 reconstruct_coupled_field 时乘入一次。
    coupled=solve_coupled_response(base,rough.response,qf;
                                   p_h=rough.params,p_f=qf.params,
                                   roughness_amplitude=1.0+0im,
                                   edge_state=zeros(ComplexF64,4),zmax=zmax)
    push!(qhf_modes,coupled.response)
    push!(qv_modes,coupled.volume_response)
    push!(qw_modes,coupled.wall_response)
    push!(source_modes,coupled.unit_source)
end

qh_at_rf=reconstruct_roughness_field(grid.kappa,qh_modes;r=rf,theta=0.0,
                                     n=n_rough,epsilon=epsilon,rf=rf,c=c,
                                     weights=grid.weights,branch=:plus)
qhf_at_rf=reconstruct_coupled_field(grid.kappa,qhf_modes;r=rf,theta=0.0,
                                    n=n_rough,epsilon=epsilon,rf=rf,c=c,
                                    weights=grid.weights,branch=:plus)

nz=length(qf.z)
function reconstructed_profile(modes,component,grid,rf,c,epsilon,n)
    nmode=length(modes)
    nz=length(VonKarmanLST._state_vector(modes[1]))÷4
    values=zeros(ComplexF64,nz)
    for iz in 1:nz
        values[iz]=reconstruct_roughness_field(grid.kappa,modes;r=rf,theta=0.0,
                    n=n,epsilon=epsilon,rf=rf,c=c,weights=grid.weights,
                    branch=:plus)[(component-1)*nz+iz]
    end
    return values
end

qh_z=reconstructed_profile(qh_modes,1,grid,rf,c,epsilon,n_rough)
qhf_z=reconstructed_profile(qhf_modes,1,grid,rf,c,epsilon,n_rough)
qv_z=reconstructed_profile(qv_modes,1,grid,rf,c,epsilon,n_rough)
qw_z=reconstructed_profile(qw_modes,1,grid,rf,c,epsilon,n_rough)
Fhf_z=reconstructed_profile(source_modes,1,grid,rf,c,epsilon,n_rough)

println("q_f residual = ",norm(qf.interior_residual))
println("q_h residual = ",maximum(norm.(getfield.(qh_modes,:interior_residual))))
println("q_hf residual = ",maximum(norm.(getfield.(qhf_modes,:interior_residual))))
println("Gaussian Hhat(0) = ",roughness_fourier_spectrum(0.0;rf=rf,c=c),
        ", kappa range = [",first(grid.kappa),", ",last(grid.kappa),"]")
mid=ceil(Int,nk/2); println("q_h wall velocity at central kappa = ",
        qh_modes[mid].q[[1,nz+1,2*nz+1]])
println("reconstructed |q_h(r_f,z=0)| = ",norm(qh_at_rf[[1,nz+1,2*nz+1]]))
println("reconstructed |q_hf(r_f,z=0)| = ",norm(qhf_at_rf[[1,nz+1,2*nz+1]]))
println("max |u_f| = ",maximum(abs,qf.q[1:nz]))
println("max reconstructed |u_h| = ",maximum(abs,qh_z))
println("max reconstructed |u_hf| = ",maximum(abs,qhf_z))
println("max reconstructed |F_hf,u| = ",maximum(abs,Fhf_z))
println("max |q_hf-(q_hf^v+q_hf^w)| = ",
        maximum(abs,reconstructed_profile(qhf_modes,1,grid,rf,c,epsilon,n_rough)-qv_z-qw_z))

# ---------------------------- 简单 SVG ----------------------------
function polyline(xs,ys,x0,y0,w,h,xmin,xmax,ymin,ymax;color="#1f77b4")
    pts=String[]
    for (x,y) in zip(xs,ys)
        xx=x0+w*(x-xmin)/(xmax-xmin)
        yy=y0+h*(1-(y-ymin)/(ymax-ymin))
        push!(pts,"$(round(xx,digits=2)),$(round(yy,digits=2))")
    end
    return "<polyline fill=\"none\" stroke=\"$color\" stroke-width=\"1.6\" points=\"$(join(pts,' '))\"/>"
end
function panel_title(x,y,s)
    return "<text x=\"$x\" y=\"$y\" font-size=\"14\" font-family=\"Arial\">$s</text>"
end
function axis_box(x,y,w,h)
    return "<rect x=\"$x\" y=\"$y\" width=\"$w\" height=\"$h\" fill=\"white\" stroke=\"#777\"/>"
end

W=1200; H=760; xpad=60; ypad=55; pw=500; ph=270
svg=String["<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$W\" height=\"$H\">",
           "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>",
           "<text x=\"40\" y=\"30\" font-size=\"20\" font-family=\"Arial\">Frozen-R coupled receptivity benchmark</text>"]
rvals=collect(range(rf-3c,rf+3c;length=161)); hvals=real.(gaussian_roughness(rvals;rf=rf,c=c,n=n_rough,branch=:real))
append!(svg,[axis_box(xpad,ypad,pw,ph),panel_title(xpad+10,ypad+22,"Gaussian roughness H_r(r)"),
             polyline(rvals,hvals,xpad,ypad,pw,ph,minimum(rvals),maximum(rvals),-0.05,1.05;color="#d62728")])
append!(svg,[axis_box(xpad+pw+80,ypad,pw,ph),panel_title(xpad+pw+90,ypad+22,"|U_h| at r=r_f"),
             polyline(qf.z,abs.(qh_z),xpad+pw+80,ypad,pw,ph,minimum(qf.z),maximum(qf.z),0,maximum(abs.(qh_z))*1.05;color="#2ca02c")])
append!(svg,[axis_box(xpad,ypad+ph+70,pw,ph),panel_title(xpad+10,ypad+ph+92,"|u_f| and |u_hf| at r=r_f"),
             polyline(qf.z,abs.(qf.q[1:nz]),xpad,ypad+ph+70,pw,ph,minimum(qf.z),maximum(qf.z),0,maximum(abs.(qf.q[1:nz]))*1.05;color="#1f77b4"),
             polyline(qf.z,abs.(qhf_z),xpad,ypad+ph+70,pw,ph,minimum(qf.z),maximum(qf.z),0,maximum(abs.(qf.q[1:nz]))*1.05;color="#9467bd")])
append!(svg,[axis_box(xpad+pw+80,ypad+ph+70,pw,ph),panel_title(xpad+pw+90,ypad+ph+92,"|F_hf,u|, |q_hf^v|, |q_hf^w|"),
             polyline(qf.z,abs.(Fhf_z),xpad+pw+80,ypad+ph+70,pw,ph,minimum(qf.z),maximum(qf.z),0,maximum(abs.(Fhf_z))*1.05;color="#ff7f0e"),
             polyline(qf.z,abs.(qv_z),xpad+pw+80,ypad+ph+70,pw,ph,minimum(qf.z),maximum(qf.z),0,maximum(abs.(qv_z))*1.05;color="#17becf"),
             polyline(qf.z,abs.(qw_z),xpad+pw+80,ypad+ph+70,pw,ph,minimum(qf.z),maximum(qf.z),0,maximum(abs.(qw_z))*1.05;color="#8c564b")])
push!(svg,"<text x=\"60\" y=\"745\" font-size=\"12\" font-family=\"Arial\">blue q_f, purple q_hf, orange F_hf, cyan volume part, brown wall part</text>")
push!(svg,"</svg>")
open(joinpath(@__DIR__,"coupled_receptivity_example.svg"),"w") do io
    write(io,join(svg,"\n"))
end
println("wrote coupled_receptivity_example.svg")

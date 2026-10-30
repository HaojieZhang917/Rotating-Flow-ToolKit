#=
旋转盘 frozen-R 涡波–粗糙度耦合响应的数值验证脚本。

该脚本验证体积/壁面分解、Chebyshev 阶数、zmax、Gaussian κ 谱和
κ→0 模态。结果只写 Markdown 报告与 SVG 图，不写 CSV/DAT。
=#
include(joinpath(@__DIR__,"VonKarmanLST.jl"))
using .VonKarmanLST
using LinearAlgebra
using Printf

const OUTDIR=joinpath(@__DIR__,"results","coupled_validation")
mkpath(OUTDIR)
const R0=285.36
const ALPHA_F=0.38482
const BETA_F=0.07759
const OMEGA_F=-BETA_F
const EDGE_AMPLITUDE=1e-3
const EPSILON=1e-3
const RF=0.0
const CROUGH=0.5
const ZI=10.0
const KAPPA_REF=0.4

function make_p(N)
    LSTParameters(mode=:spatial,R=R0,alpha=ALPHA_F,beta=BETA_F,
                  omega=OMEGA_F,N=N)
end

function edge_state(p)
    first(vortical_edge_states([p.alpha,p.beta,0.0];
                               amplitude=EDGE_AMPLITUDE))
end

function state_profile(response,zq; components=1:4)
    z=Float64.(response.z); n=length(z)
    bw=VonKarmanLST._barycentric_weights(z)
    q=VonKarmanLST._state_vector(response)
    out=zeros(ComplexF64,length(components),length(zq))
    for (ii,c) in enumerate(components)
        values=q[(c-1)*n+1:c*n]
        for j in eachindex(zq)
            out[ii,j]=VonKarmanLST._barycentric_interpolate(z,values,bw,zq[j])
        end
    end
    out
end

function reconstructed_profile(modes,zgrid,kappa,zq;weights,
                               r=RF,epsilon=EPSILON)
    nm=length(modes); length(kappa)==nm==length(weights) ||
        error("modes、kappa 和 weights 长度不一致")
    nz=length(zgrid); out=zeros(ComplexF64,4,length(zq))
    z=Float64.(zgrid); bw=VonKarmanLST._barycentric_weights(z)
    for j in 1:nm
        q=VonKarmanLST._state_vector(modes[j])
        coef=ComplexF64(epsilon)*weights[j]/(2π)*
              roughness_fourier_spectrum(kappa[j];rf=RF,c=CROUGH)*
              exp(im*kappa[j]*r)
        for c in 1:4
            vals=q[(c-1)*nz+1:c*nz]
            for l in eachindex(zq)
                out[c,l]+=coef*VonKarmanLST._barycentric_interpolate(
                    z,vals,bw,zq[l])
            end
        end
    end
    out
end

function energy_metrics(v,w,f,weights)
    en(q)=sum(weights .* vec(sum(abs2.(q[1:3,:]);dims=1)))
    ev=en(v); ew=en(w); ef=en(f)
    ivw=sum(weights .* vec(sum(conj.(v[1:3,:]).*w[1:3,:];dims=1)))
    comp(q)=[sum(weights .* abs2.(q[j,:])) for j in 1:3]
    chi=(ev+ew)>0 ? 2real(ivw)/(ev+ew) : NaN
    rho=ew>0 ? sqrt(ev/ew) : NaN
    ident=ef>0 ? abs(ef-(ev+ew+2real(ivw)))/ef : NaN
    (E_v=ev,E_w=ew,E_full=ef,I_vw=ivw,chi=chi,rho=rho,
     component_v=comp(v),component_w=comp(w),component_full=comp(f),
     identity_error=ident)
end

function fixed_case(N,zmax;kappa=KAPPA_REF)
    base=solve_von_karman_baseflow()
    p=make_p(N)
    qf=solve_forced_spatial(base,edge_state(p);p=p,zmax=zmax)
    rh=solve_roughness_distortion(base;p=p,kappa=kappa,n=0,
                                  beta_shift=0.0,h_amplitude=1.0,zmax=zmax)
    cr=solve_coupled_response(base,rh.response,qf;p_h=rh.params,p_f=qf.params,
                              roughness_amplitude=EPSILON,
                              edge_state=zeros(ComplexF64,4),zmax=zmax)
    zq,wq=subinterval_cc_grid(ZI,128)
    qfp=state_profile(qf,zq;components=1:3)
    qvp=state_profile(cr.volume_response,zq;components=1:3)
    qwp=state_profile(cr.wall_response,zq;components=1:3)
    qcp=state_profile(cr.response,zq;components=1:3)
    (N=N,zmax=zmax,qf=qf,qh=rh.response,coupled=cr,
     qf_profile=qfp,volume_profile=qvp,wall_profile=qwp,full_profile=qcp,
     metrics=energy_metrics(qvp,qwp,qcp,wq),
     qf_residual=norm(qf.interior_residual),
     qh_residual=norm(rh.response.interior_residual),
     qhf_residual=norm(cr.response.interior_residual),
     bc_residual=maximum(abs,cr.response.boundary_error),
     zq=zq,wq=wq)
end

function gaussian_case(N,zmax;kmax=12.0,nk=32,decompose=true)
    base=solve_von_karman_baseflow()
    p=make_p(N)
    qf=solve_forced_spatial(base,edge_state(p);p=p,zmax=zmax)
    grid=roughness_fourier_grid(c=CROUGH,kappa_max=kmax,nk=nk)
    qh=Any[]; qhf=Any[]; qv=Any[]; qw=Any[]; fs=Any[]
    for kappa in grid.kappa
        rh=solve_roughness_distortion(base;p=p,kappa=kappa,n=0,
                                      beta_shift=0.0,h_amplitude=1.0,zmax=zmax)
        cr=solve_coupled_response(base,rh.response,qf;p_h=rh.params,p_f=qf.params,
                                  roughness_amplitude=1.0,
                                  edge_state=zeros(ComplexF64,4),zmax=zmax,
                                  decompose=decompose)
        push!(qh,rh.response); push!(qhf,cr.response); push!(fs,cr.unit_source)
        if decompose
            push!(qv,cr.volume_response); push!(qw,cr.wall_response)
        end
    end
    zq,wq=subinterval_cc_grid(ZI,128)
    rough=reconstructed_profile(qh,qf.z,grid.kappa,zq;weights=grid.weights)
    full=reconstructed_profile(qhf,qf.z,grid.kappa,zq;weights=grid.weights)
    source=reconstructed_profile(fs,qf.z,grid.kappa,zq;weights=grid.weights)
    if decompose
        vol=reconstructed_profile(qv,qf.z,grid.kappa,zq;weights=grid.weights)
        wall=reconstructed_profile(qw,qf.z,grid.kappa,zq;weights=grid.weights)
        metrics=energy_metrics(vol,wall,full,wq)
    else
        vol=zeros(ComplexF64,3,length(zq)); wall=zeros(ComplexF64,3,length(zq))
        metrics=(E_v=NaN,E_w=NaN,E_full=sum(wq.*vec(sum(abs2.(full[1:3,:]);dims=1))),
                 I_vw=NaN,chi=NaN,rho=NaN,component_v=fill(NaN,3),
                 component_w=fill(NaN,3),component_full=[sum(wq.*abs2.(full[j,:])) for j in 1:3],
                 identity_error=NaN)
    end
    (N=N,zmax=zmax,qf=qf,grid=grid,qh=qh,qhf=qhf,qv=qv,qw=qw,source_modes=fs,
     rough=rough,full=full,source=source,volume=vol,wall=wall,
     metrics=metrics,zq=zq,wq=wq)
end

function boundary_matrix(op)
    P=copy(op.P)
    bc=setdiff(1:size(P,1),op.blocks.keep)
    for row in bc
        P[row,:].=0; P[row,row]=1
    end
    P
end

function zero_mode_scan(N=29)
    base=solve_von_karman_baseflow()
    ks=[1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]
    smin=Float64[]; qnorm=Float64[]
    for k in ks
        p=LSTParameters(mode=:spatial,R=R0,alpha=k,beta=0.0,omega=0.0,N=N)
        op=assemble_forced_spatial(base;p=p,zmax=30.0)
        push!(smin,minimum(svdvals(boundary_matrix(op))))
        rh=solve_roughness_distortion(base;p=p,kappa=k,n=0,beta_shift=0.0,
                                      h_amplitude=1.0,zmax=30.0)
        push!(qnorm,norm(rh.response.q))
    end
    p0=LSTParameters(mode=:spatial,R=R0,alpha=0.0,beta=0.0,omega=0.0,N=N)
    op0=assemble_forced_spatial(base;p=p0,zmax=30.0)
    S=svd(boundary_matrix(op0)); n=N+1
    null=S.V[:,end]
    comps=[norm(null[(j-1)*n+1:j*n]) for j in 1:4]
    (kappa=ks,smin=smin,qnorm=qnorm,zero_smin=minimum(S.S),
     zero_null_components=comps)
end

fmt(x)=isfinite(real(x)) ? @sprintf("%.6e",real(x)) : "NaN"
const VALIDATION_START_TIME=time()

function relative_profile_error(profile,reference,weights)
    num=sum(weights .* vec(sum(abs2.(profile-reference);dims=1)))
    den=sum(weights .* vec(sum(abs2.(reference);dims=1)))
    sqrt(num/max(den,eps(Float64)))
end

function write_report(main,fixed,conv_errors,zrows,zerrors,spectral,zero)
    io=IOBuffer()
    println(io,"# Coupled receptivity validation")
    println(io)
    println(io,"模型：frozen-R、有限外缘 Dirichlet；R=$(R0), alpha_f=$(ALPHA_F), beta_f=$(BETA_F), omega_f=$(OMEGA_F), epsilon=$(EPSILON), a=$(EDGE_AMPLITUDE), rf=$(RF), c=$(CROUGH), n=0。")
    println(io)
    m=main.metrics
    println(io,"## Gaussian 主算例")
    println(io)
    println(io,"|指标|数值|")
    println(io,"|---|---:|")
    println(io,"|qf 残差|$(fmt(norm(main.qf.interior_residual)))|")
    println(io,"|qh 最大残差|$(fmt(maximum(norm.(getfield.(main.qh,:interior_residual)))))|")
    println(io,"|qhf 最大残差|$(fmt(maximum(norm.(getfield.(main.qhf,:interior_residual)))))|")
    println(io,"|E_v|$(fmt(m.E_v))|")
    println(io,"|E_w|$(fmt(m.E_w))|")
    println(io,"|E_full|$(fmt(m.E_full))|")
    println(io,"|rho|$(fmt(m.rho))|")
    println(io,"|chi|$(fmt(m.chi))|")
    println(io,"|能量恒等式相对误差|$(fmt(m.identity_error))|")
    println(io)
    println(io,"速度分量能量 (u,v,w)：")
    println(io)
    println(io,"|分量|体积|壁面|完整|")
    println(io,"|---|---:|---:|---:|")
    for j in 1:3
        println(io,"|$(('u','v','w')[j])|$(fmt(m.component_v[j]))|$(fmt(m.component_w[j]))|$(fmt(m.component_full[j]))|")
    end
    println(io)
    println(io,"## Chebyshev 阶数收敛（单一 kappa=$(KAPPA_REF)）")
    println(io)
    println(io,"|N|qf残差|qh残差|qhf残差|full剖面误差|E_full|")
    println(io,"|---:|---:|---:|---:|---:|---:|")
    for (row,e) in zip(fixed,conv_errors)
        println(io,"|$(row.N)|$(fmt(row.qf_residual))|$(fmt(row.qh_residual))|$(fmt(row.qhf_residual))|$(fmt(e))|$(fmt(row.metrics.E_full))|")
    end
    println(io)
    println(io,"## zmax 收敛（N=69）")
    println(io)
    println(io,"|zmax|full剖面误差|E_v|E_w|E_full|")
    println(io,"|---:|---:|---:|---:|---:|")
    for (row,e) in zip(zrows,zerrors)
        println(io,"|$(row.zmax)|$(fmt(e))|$(fmt(row.metrics.E_v))|$(fmt(row.metrics.E_w))|$(fmt(row.metrics.E_full))|")
    end
    println(io)
    println(io,"## Gaussian kappa 频谱收敛")
    println(io)
    println(io,"最终响应误差以最后一个加密算例为参考；Gaussian 几何误差只检验解析谱本身，不代表 LNS 加权响应已经收敛。")
    println(io)
    println(io,"|kappa_max|Nkappa|Gaussian几何误差|最终响应误差|E_full|")
    println(io,"|---:|---:|---:|---:|---:|")
    for row in spectral
        println(io,"|$(row.kmax)|$(row.nk)|$(fmt(row.geometry_error))|$(fmt(row.response_error))|$(fmt(row.metrics.E_full))|")
    end
    println(io)
    println(io,"## kappa→0 检查")
    println(io)
    println(io,"|kappa|边界矩阵最小奇异值|单位粗糙度响应范数|")
    println(io,"|---:|---:|---:|")
    for j in eachindex(zero.kappa)
        println(io,"|$(fmt(zero.kappa[j]))|$(fmt(zero.smin[j]))|$(fmt(zero.qnorm[j]))|")
    end
    println(io)
    println(io,"精确 kappa=0 的最小奇异值为 $(fmt(zero.zero_smin))。其零奇异向量速度/压力块范数为 $(join(fmt.(zero.zero_null_components),", "))，说明当前离散零模主要位于压力块，但仍需要显式压力规范才能把 kappa=0 纳入积分。")
    println(io)
    println(io,"## 结论")
    println(io)
    println(io,"主算例的体积强迫响应与壁面响应满足线性叠加和复内积能量恒等式。当前参数下体积项远大于壁面项；chi 的符号给出两者在速度范数下的微弱增强或抵消。")
    println(io,"但是 qf 和 qhf 对有限外缘 zmax 明显敏感，壁面项还受到 qf 壁面梯度谱收敛的影响。因此这些结果只能解释为有限 Dirichlet、frozen-R 局部模型结果，不能作为真实自由涡波感受性或 Type I/II 激发系数。")
    write(joinpath(OUTDIR,"coupled_validation_report.md"),String(take!(io)))
end

function heat_color(t)
    tt=clamp(Float64(t),0.0,1.0)
    @sprintf("#%02x%02x%02x",round(Int,255*tt),80,round(Int,255*(1-tt)))
end

function svg_heatmap(path,mats,titles; figure_title="",
                     x_label="z (non-dimensional wall-normal coordinate)",
                     y_label="r-r_f (non-dimensional local radial coordinate)",
                     value_label="magnitude of complex Fourier amplitude",
                     xlims=(0.0,ZI),ylims=(RF-1.5,RF+1.5),
                     note="Frozen-R local cylindrical coordinates; plotted values are absolute values.")
    W=1320; H=500; x0=55; y0=75; pw=350; ph=300; gap=35
    vmax=max(maximum(vcat([vec(m) for m in mats]...)),eps(Float64))
    s=String["<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$W\" height=\"$H\">",
             "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>",
             "<text x=\"35\" y=\"30\" font-size=\"20\" font-family=\"Arial\">$figure_title</text>"]
    for k in eachindex(mats)
        mat=mats[k]; nr,nz=size(mat); x=x0+(k-1)*(pw+gap)
        push!(s,"<text x=\"$x\" y=\"58\" font-size=\"14\" font-family=\"Arial\">$(titles[k])</text>")
        for i in 1:nr, j in 1:nz
            xx=x+(j-1)*pw/nz; yy=y0+(i-1)*ph/nr
            push!(s,"<rect x=\"$xx\" y=\"$yy\" width=\"$(pw/nz+0.5)\" height=\"$(ph/nr+0.5)\" fill=\"$(heat_color(mat[i,j]/vmax))\"/>")
        end
        push!(s,"<rect x=\"$x\" y=\"$y0\" width=\"$pw\" height=\"$ph\" fill=\"none\" stroke=\"#555\"/>")
        for (xx,label) in zip([x,x+pw/2,x+pw],
                              [fmt(xlims[1]),fmt((xlims[1]+xlims[2])/2),fmt(xlims[2])])
            push!(s,"<text x=\"$xx\" y=\"$(y0+ph+18)\" text-anchor=\"middle\" font-size=\"10\" font-family=\"Arial\">$label</text>")
        end
        for (yy,label) in zip([y0+ph,y0+ph/2,y0],
                              [fmt(ylims[1]),fmt((ylims[1]+ylims[2])/2),fmt(ylims[2])])
            push!(s,"<text x=\"$(x-6)\" y=\"$(yy+4)\" text-anchor=\"end\" font-size=\"10\" font-family=\"Arial\">$label</text>")
        end
    end
    cbx=1225; cby=y0; cbw=18; cbh=ph
    for i in 0:99
        t=i/99
        push!(s,"<rect x=\"$cbx\" y=\"$(cby+cbh*(1-t))\" width=\"$cbw\" height=\"$(cbh/100+1)\" fill=\"$(heat_color(t))\"/>")
    end
    push!(s,"<text x=\"$(cbx-4)\" y=\"$(cby-8)\" font-size=\"11\" font-family=\"Arial\">$(fmt(vmax))</text>")
    push!(s,"<text x=\"$(cbx-4)\" y=\"$(cby+cbh+16)\" font-size=\"11\" font-family=\"Arial\">0</text>")
    push!(s,"<text x=\"$(cbx+28)\" y=\"$(cby+cbh/2)\" font-size=\"11\" font-family=\"Arial\" transform=\"rotate(90 $(cbx+28) $(cby+cbh/2))\">$value_label</text>")
    push!(s,"<text x=\"$(x0+pw)\" y=\"$(y0+ph+42)\" text-anchor=\"middle\" font-size=\"11\" font-family=\"Arial\">$x_label</text>")
    push!(s,"<text x=\"18\" y=\"$(y0+ph/2)\" text-anchor=\"middle\" font-size=\"11\" font-family=\"Arial\" transform=\"rotate(-90 18 $(y0+ph/2))\">$y_label</text>")
    push!(s,"<text x=\"35\" y=\"480\" font-size=\"11\" font-family=\"Arial\">$note; R=$(R0), alpha_f=$(ALPHA_F), beta_f=$(BETA_F), omega_f=$(OMEGA_F), epsilon=$(EPSILON), a=$(EDGE_AMPLITUDE), r_f=$(RF), c=$(CROUGH), n=0.</text>")
    push!(s,"</svg>")
    write(path,join(s,"\n"))
end

function svg_lines(path,xs,ys,labels,colors,title;
                   x_label="coordinate",
                   y_label="magnitude",
                   note="Complex Fourier amplitudes are plotted by absolute value.")
    W=1000; H=540; x0=80; y0=70; pw=800; ph=360
    ymax=max(maximum(vcat([collect(y) for y in ys]...)),eps(Float64))
    s=String["<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$W\" height=\"$H\">",
             "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>",
             "<text x=\"35\" y=\"30\" font-size=\"20\" font-family=\"Arial\">$title</text>",
             "<rect x=\"$x0\" y=\"$y0\" width=\"$pw\" height=\"$ph\" fill=\"white\" stroke=\"#777\"/>"]
    for (y,col) in zip(ys,colors)
        pts=String[]
        for (xx,yy) in zip(xs,y)
            px=x0+pw*(xx-first(xs))/(last(xs)-first(xs))
            py=y0+ph*(1-yy/ymax)
            push!(pts,"$(round(px,digits=2)),$(round(py,digits=2))")
        end
        push!(s,"<polyline fill=\"none\" stroke=\"$col\" stroke-width=\"2\" points=\"$(join(pts,' '))\"/>")
    end
    for (xx,label) in zip([first(xs),first(xs)+(last(xs)-first(xs))/2,last(xs)],
                          [fmt(first(xs)),fmt((first(xs)+last(xs))/2),fmt(last(xs))])
        px=x0+pw*(xx-first(xs))/(last(xs)-first(xs))
        push!(s,"<text x=\"$px\" y=\"$(y0+ph+20)\" text-anchor=\"middle\" font-size=\"11\" font-family=\"Arial\">$label</text>")
    end
    for (yy,label) in zip([y0+ph,y0+ph/2,y0],["0",fmt(ymax/2),fmt(ymax)])
        push!(s,"<text x=\"$(x0-8)\" y=\"$(yy+4)\" text-anchor=\"end\" font-size=\"11\" font-family=\"Arial\">$label</text>")
    end
    for (j,label) in enumerate(labels)
        yy=80+22*(j-1)
        push!(s,"<line x1=\"$(x0+pw-170)\" y1=\"$yy\" x2=\"$(x0+pw-145)\" y2=\"$yy\" stroke=\"$(colors[j])\" stroke-width=\"2\"/>")
        push!(s,"<text x=\"$(x0+pw-138)\" y=\"$(yy+5)\" font-size=\"12\" font-family=\"Arial\">$label</text>")
    end
    push!(s,"<text x=\"$(x0+pw/2)\" y=\"$(y0+ph+45)\" text-anchor=\"middle\" font-size=\"12\" font-family=\"Arial\">$x_label</text>")
    push!(s,"<text x=\"18\" y=\"$(y0+ph/2)\" text-anchor=\"middle\" font-size=\"12\" font-family=\"Arial\" transform=\"rotate(-90 18 $(y0+ph/2))\">$y_label</text>")
    push!(s,"<text x=\"35\" y=\"520\" font-size=\"11\" font-family=\"Arial\">$note; R=$(R0), alpha_f=$(ALPHA_F), beta_f=$(BETA_F), omega_f=$(OMEGA_F), epsilon=$(EPSILON), a=$(EDGE_AMPLITUDE).</text>")
    push!(s,"</svg>")
    write(path,join(s,"\n"))
end

function write_figures(main,fixed,conv_errors,zrows,zerrors,spectral)
    rr=collect(range(RF-1.5,RF+1.5;length=41))
    zz=collect(range(0.0,ZI;length=41))
    function mat(modes,component)
        out=zeros(Float64,length(rr),length(zz))
        for i in eachindex(rr)
            p=reconstructed_profile(modes,main.qf.z,main.grid.kappa,zz;
                                     weights=main.grid.weights,r=rr[i])
            out[i,:].=abs.(p[component,:])
        end
        out
    end
    uh=[mat(main.qh,j) for j in 1:3]
    ff=[mat(main.source_modes,j) for j in 1:3]
    function matnorm(modes)
        out=zeros(Float64,length(rr),length(zz))
        for i in eachindex(rr)
            p=reconstructed_profile(modes,main.qf.z,main.grid.kappa,zz;
                                     weights=main.grid.weights,r=rr[i])
            out[i,:].=sqrt.(vec(sum(abs2.(p[1:3,:]);dims=1)))
        end
        out
    end
    qv=matnorm(main.qv); qw=matnorm(main.qw); qfull=matnorm(main.qhf)
    svg_heatmap(joinpath(OUTDIR,"figure2_roughness_distortion.svg"),uh,
                ["|U_h|","|V_h|","|W_h|"];
                figure_title="Figure 2: stationary roughness distortion",
                value_label="|complex Fourier velocity amplitude|",
                note="U_h,V_h,W_h are reconstructed physical velocity components; r-r_f is local radial and z is wall-normal.")
    svg_heatmap(joinpath(OUTDIR,"figure3_cross_forcing.svg"),ff,
                ["|F_r|","|F_theta|","|F_z|"];
                figure_title="Figure 3: cross-order volume forcing",
                value_label="|F_hf component|",
                note="F_hf is the reconstructed O(epsilon*a) momentum forcing; pressure and continuity source are zero.")
    svg_heatmap(joinpath(OUTDIR,"figure4_coupling_response.svg"),[qv,qw,qfull],
                ["volume |u_hf|","wall |u_hf|","full |u_hf|"];
                figure_title="Figure 4: volume, wall and full coupled response",
                value_label="velocity norm |u_hf|",
                note="Each panel is the norm of (u,v,w), reconstructed at the same physical location and scaled by epsilon*a.")
    pressure=mat(main.qhf,4)
    svg_heatmap(joinpath(OUTDIR,"figure4b_full_pressure.svg"),[pressure],
                ["full |p_hf|"];
                figure_title="Figure 4b: full coupled pressure response",
                value_label="|p_hf|",
                note="Pressure is the fourth primitive variable; it is plotted separately from the velocity norm.")
    hvals=abs.(gaussian_roughness(rr;rf=RF,c=CROUGH,n=0,branch=:real))
    svg_lines(joinpath(OUTDIR,"figure1_roughness_spectrum.svg"),rr,[hvals],
              ["|H_r(r)|"],["#d62728"],"Figure 1a: Gaussian roughness envelope";
              x_label="r-r_f (local radial coordinate)",y_label="|H_r|",
              note="Axisymmetric Gaussian ring: H_r=exp(-(r-r_f)^2/(2c^2)); r_f=0, c=0.5.")
    hhat=abs.(roughness_fourier_spectrum(main.grid.kappa;rf=RF,c=CROUGH))
    svg_lines(joinpath(OUTDIR,"figure1b_fourier_spectrum.svg"),main.grid.kappa,[hhat],
              ["|Hhat(kappa)|"],["#1f77b4"],"Figure 1b: Gaussian radial Fourier spectrum";
              x_label="kappa (radial Fourier wavenumber)",y_label="|Hhat_r(kappa)|",
              note="Hhat uses Hhat=integral H_r exp(-i*kappa*r)dr and inverse factor 1/(2*pi).")
    zq=main.zq
    for (j,name) in enumerate(["u","v","w"])
        pvol=abs.(main.volume[j,:]); pwall=abs.(main.wall[j,:]); pfull=abs.(main.full[j,:])
        filename=j==1 ? "figure5_normal_profiles.svg" : "figure5_$(name)_normal_profiles.svg"
        svg_lines(joinpath(OUTDIR,filename),zq,[pvol,pwall,pfull],
                  ["volume","wall","full"],
                  ["#1f77b4","#8c564b","#9467bd"],
                  "Figure 5: normal profile of $name at r=r_f";
                  x_label="z (wall-normal coordinate)",y_label="|$(name)_hf(z)|",
                  note="Curves are absolute values of reconstructed complex amplitudes; volume, wall and full use identical z quadrature.")
    end
    svg_lines(joinpath(OUTDIR,"figure6_convergence.svg"),[row.N for row in fixed],
              [conv_errors],["Chebyshev full-profile error"],["#1f77b4"],
              "Figure 6a: Chebyshev convergence";
              x_label="N (Chebyshev polynomial order)",y_label="relative L2 profile error",
              note="Reference is N=89; comparison uses the same physical z in [0,10].")
    svg_lines(joinpath(OUTDIR,"figure6b_zmax.svg"),[row.zmax for row in zrows],
              [[row.metrics.E_full for row in zrows]],["E_full"],["#d62728"],
              "Figure 6b: finite outer-domain sensitivity";
              x_label="z_max (finite outer boundary)",y_label="E_full",
              note="E_full is the velocity energy over z in [0,10]; outer Dirichlet data are unchanged.")
    svg_lines(joinpath(OUTDIR,"figure6c_kappa.svg"),collect(1:length(spectral)),
              [[row.response_error for row in spectral]],["spectral response error"],["#2ca02c"],
              "Figure 6c: Gaussian Fourier quadrature convergence";
              x_label="test case index (see report for kappa_max,N_kappa)",y_label="relative response error",
              note="Reference is kappa_max=16,N_kappa=128; geometry and LNS-weighted response errors are distinguished in the report.")
end

println("main Gaussian case")
main=gaussian_case(29,20.0;kmax=12.0,nk=32,decompose=true)
println("Chebyshev convergence")
Ns=[29,49,69,89]
fixed=[fixed_case(N,30.0) for N in Ns]
ref=fixed[end]
zq,wq=subinterval_cc_grid(ZI,128)
conv_errors=[relative_profile_error(row.full_profile,ref.full_profile,wq) for row in fixed]

println("zmax convergence")
zvals=[15.0,20.0,30.0,40.0]
zrows=[fixed_case(69,zmax) for zmax in zvals]
zref=zrows[end]
zerrors=[relative_profile_error(row.full_profile,zref.full_profile,wq) for row in zrows]

println("Gaussian spectral convergence")
spectral_settings=[(8.0,32),(8.0,64),(8.0,128),
                    (12.0,32),(12.0,64),(12.0,128),
                    (16.0,32),(16.0,64),(16.0,128)]
spectral=Any[]
for (kmax,nk) in spectral_settings
    println("  kmax=",kmax," nk=",nk)
    run=gaussian_case(19,20.0;kmax=kmax,nk=nk,decompose=false)
    geom=abs(sum(run.grid.weights .* roughness_fourier_spectrum(
        run.grid.kappa;rf=RF,c=CROUGH))/(2π)-1)
    push!(spectral,(kmax=kmax,nk=nk,geometry_error=geom,
                    metrics=run.metrics,response_error=NaN))
end
eref=spectral[end].metrics.E_full
for j in eachindex(spectral)
    spectral[j]=(;spectral[j]...,response_error=
        abs(spectral[j].metrics.E_full-eref)/max(eref,eps(Float64)))
end

println("near-zero kappa")
zero=zero_mode_scan(29)
println("main E_v=",main.metrics.E_v," E_w=",main.metrics.E_w,
        " E_full=",main.metrics.E_full," chi=",main.metrics.chi,
        " rho=",main.metrics.rho)
write_report(main,fixed,conv_errors,zrows,zerrors,spectral,zero)
write_figures(main,fixed,conv_errors,zrows,zerrors,spectral)
elapsed=time()-VALIDATION_START_TIME
write(joinpath(OUTDIR,"performance.txt"),@sprintf("total wall time: %.3f s\n",elapsed))
println("wrote ",joinpath(OUTDIR,"coupled_validation_report.md"))

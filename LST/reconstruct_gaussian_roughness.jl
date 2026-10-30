include(joinpath(@__DIR__,"VonKarmanLST.jl")); using .VonKarmanLST
include(joinpath(@__DIR__,"solve_six_state_lns.jl")); using .SixStateLNS
using LinearAlgebra
using Serialization

const OUTDIR=joinpath(@__DIR__,"results","roughness_reconstruction"); mkpath(OUTDIR)
const R0=285.36; const ALPHA=.38482; const BETA=.07759; const OMEGA=-BETA
const AINC=1e-3+0im; const EPS=1e-3; const RF=0.0; const CROUGH=1.0
const ZI=10.0; const NMAIN=29; const ZMAX=20.0; const KMAX=6.0; const NK=32

function incident_state(base,p,z)
    modes=farfield_six_state_modes(base;p=p); c=[m for m in modes if imag(m.kz)<-1e-8]
    lm=minimum(abs(imag(m.kz)) for m in c); c=[m for m in c if abs(abs(imag(m.kz))-lm)<1e-8]
    m=c[argmax([real(x.kz) for x in c])]; y=m.state/norm(m.state[1:3])
    AINC.*y.*exp(im*m.kz*(z-10.0))
end

function interp_profile(resp,zq)
    z=Float64.(resp.z); n=length(z); bw=VonKarmanLST._barycentric_weights(z)
    out=zeros(ComplexF64,4,length(zq))
    for c in 1:4
        vals=ComplexF64.(resp.q[(c-1)*n+1:c*n])
        for j in eachindex(zq); out[c,j]=VonKarmanLST._barycentric_interpolate(z,vals,bw,zq[j]); end
    end
    out
end

function write_csv(path,header,rows)
    io=IOBuffer(); println(io,join(header,",")); for row in rows; println(io,join(string.(row),",")); end
    write(path,String(take!(io)))
end

function solve_spectral_case(; N=NMAIN,zmax=ZMAX,kmax=KMAX,nk=NK,save=false)
    base=solve_von_karman_baseflow(); p=LSTParameters(mode=:spatial,R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,N=N,nev=8,spatial_shift=.38+0im)
    yi=incident_state(base,p,zmax); first=solve_six_state_lns(base;p=p,zmax=zmax,boundary=:modal,incident_state=yi); qf=as_four_state_response(first)
    grid=roughness_fourier_grid(c=CROUGH,kappa_max=kmax,nk=nk)
    zq,wz=subinterval_cc_grid(ZI,128); profiles=Matrix{ComplexF64}[]; unit_energy=Float64[]; modes=Any[]
    for (j,k) in enumerate(grid.kappa)
        println("spectral mode ",j,"/",length(grid.kappa)," kappa=",k)
        qh=solve_roughness_distortion(base;p=p,kappa=k,n=0,h_amplitude=1.0,zmax=zmax,edge_state=zeros(ComplexF64,4))
        co=solve_six_state_coupled_response(base,qh.response,qf;roughness_amplitude=1.0,zmax=zmax,boundary=:modal)
        prof=interp_profile(as_four_state_response(co.response),zq); push!(profiles,prof); push!(unit_energy,sum(wz.*vec(sum(abs2.(prof[1:3,:]);dims=1)))); push!(modes,(qh=qh,co=co))
    end
    hhat=roughness_fourier_spectrum(grid.kappa;rf=RF,c=CROUGH)
    coef=EPS .* grid.weights .* hhat ./ (2π)
    Lr=2π/(grid.kappa[2]-grid.kappa[1]); nr=nk; rgrid=collect(range(-Lr/2,step=Lr/nr,length=nr)); qphys=zeros(ComplexF64,4,nr,length(zq));
    for ir in eachindex(rgrid), j in eachindex(grid.kappa); qphys[:,ir,:] .+= (coef[j]*exp(im*(ALPHA+grid.kappa[j])*rgrid[ir])).*profiles[j]; end
    local_energy=[sum(wz.*vec(sum(abs2.(qphys[1:3,ir,:]);dims=1))) for ir in 1:nr]
    dr=Lr/nr; Ephys=dr*sum(local_energy); Eparse=Lr*sum(abs2.(coef).*unit_energy)
    wall_err=Float64[]; geometry_err=Float64[]
    for ir in eachindex(rgrid)
        h=gaussian_roughness(rgrid[ir];rf=RF,c=CROUGH,epsilon=EPS,branch=:real)
        hspec=sum(coef.*exp.(im.*grid.kappa.*rgrid[ir])); push!(geometry_err,abs(hspec-h))
        # q_f 的径向 Fourier 因子 e^{i alpha_f r} 也必须保留。
        wall_ref=-hspec.*ComplexF64.(qf.wall_shear).*exp(im*ALPHA*rgrid[ir])
        push!(wall_err,norm(qphys[1:3,ir,1]-wall_ref)/max(norm(wall_ref),eps(Float64)))
    end
    single_gain=unit_energy./abs(AINC)^2; weighted=abs2.(hhat).*single_gain
    if save
        rows=[]
        for ir in eachindex(rgrid), iz in eachindex(zq)
            push!(rows,(rgrid[ir],zq[iz],real(qphys[1,ir,iz]),imag(qphys[1,ir,iz]),real(qphys[2,ir,iz]),imag(qphys[2,ir,iz]),real(qphys[3,ir,iz]),imag(qphys[3,ir,iz]),real(qphys[4,ir,iz]),imag(qphys[4,ir,iz]),local_energy[ir]))
        end
        write_csv(joinpath(OUTDIR,"reconstructed_response.csv"),["r","z","u_re","u_im","v_re","v_im","w_re","w_im","p_re","p_im","E_local"],rows)
        write_csv(joinpath(OUTDIR,"local_energy.csv"),["r","E_hf"],[(rgrid[i],local_energy[i]) for i in eachindex(rgrid)])
        write_csv(joinpath(OUTDIR,"roughness_spectrum.csv"),["kappa","Hhat_re","Hhat_im","abs_Hhat","G_single","abs_Hhat_sq_G_single"],[(grid.kappa[i],real(hhat[i]),imag(hhat[i]),abs(hhat[i]),single_gain[i],weighted[i]) for i in eachindex(grid.kappa)])
        serialize(joinpath(OUTDIR,"reconstructed_response.bin"),(r=rgrid,z=zq,q=qphys,local_energy=local_energy,kappa=grid.kappa,weights=grid.weights,coefficients=coef,E_physical=Ephys,E_parseval=Eparse,wall_error=wall_err,geometry_error=geometry_err,params=(R=R0,alpha=ALPHA,beta=BETA,omega=OMEGA,epsilon=EPS,rf=RF,c=CROUGH,N=N,zmax=zmax)))
    end
    return (grid=grid,r=rgrid,z=zq,q=qphys,profiles=profiles,hhat=hhat,coef=coef,unit_energy=unit_energy,single_gain=single_gain,weighted=weighted,local_energy=local_energy,Ephys=Ephys,Eparse=Eparse,wall_error=wall_err,qf=qf,modes=modes,Lr=Lr)
end

function gaussian_validation()
    rr=collect(range(-3CROUGH,3CROUGH;length=121)); g=roughness_fourier_grid(c=CROUGH,kappa_max=KMAX,nk=256); hhat=roughness_fourier_spectrum(g.kappa;rf=RF,c=CROUGH);
    hrec=ComplexF64[]; for r in rr; push!(hrec,sum(g.weights.*hhat.*exp.(im.*g.kappa.*r))/(2π)); end
    hex=gaussian_roughness(rr;rf=RF,c=CROUGH,epsilon=1.0,branch=:real); err=abs.(hrec-hex)
    write_csv(joinpath(OUTDIR,"fourier_transform_validation.csv"),["r","h_exact","h_reconstructed","absolute_error"],[(rr[i],real(hex[i]),real(hrec[i]),err[i]) for i in eachindex(rr)])
    return (max_error=maximum(err),peak_exact=maximum(real.(hex)),peak_reconstructed=maximum(real.(hrec)),peak_r=rr[argmax(real.(hrec))])
end

println("validating analytic Gaussian Fourier transform")
fv=gaussian_validation(); println(fv)
println("solving main six-state Modal spectral reconstruction")
main=solve_spectral_case(save=true)
    println("E_physical=",main.Ephys," E_parseval=",main.Eparse," relative=",abs(main.Ephys-main.Eparse)/main.Eparse," wall_error=",maximum(main.wall_error)," geometry_error=",maximum(main.geometry_error))
conv=[]
for (kmax,nk) in [(4.0,16),(4.0,32),(6.0,16),(6.0,32),(8.0,32),(6.0,48)]
    println("convergence kmax=",kmax," nk=",nk)
    x=solve_spectral_case(N=19,zmax=20.0,kmax=kmax,nk=nk,save=false)
    push!(conv,(kmax=kmax,nk=nk,geometry_error=abs(sum(x.grid.weights.*x.hhat)/(2π)-1),E_physical=x.Ephys,E_parseval=x.Eparse))
end
eref=main.Ephys; write_csv(joinpath(OUTDIR,"reconstruction_convergence.csv"),["kappa_max","nk","geometry_error","E_physical","E_parseval","relative_to_main"],[(x.kmax,x.nk,x.geometry_error,x.E_physical,x.E_parseval,abs(x.E_physical-eref)/eref) for x in conv])
open(joinpath(OUTDIR,"roughness_reconstruction_report.md"),"w") do io
    println(io,"# Gaussian 局部粗糙元 Fourier 重构验证")
    println(io,"\n参数：`R=$(R0)`, `alpha_f=$(ALPHA)`, `beta_f=$(BETA)`, `omega_f=$(OMEGA)`, `A_inc=$(AINC)`, `epsilon=$(EPS)`, `r_f=$(RF)`, `c=$(CROUGH)`, `N=$(NMAIN)`, `zmax=$(ZMAX)`。每个单 Fourier 响应均使用六维 Modal 二阶求解器，目标波数为 `alpha_hf=alpha_f+kappa`。")
    println(io,"\n## Fourier 变换验证")
    println(io,"- 解析 Gaussian 峰值：`1`；数值逆变换峰值：`$(fv.peak_reconstructed)`，峰值位置：`$(fv.peak_r)`。")
    println(io,"- 最大绝对重构误差：`$(fv.max_error)`。逆变换使用 `1/(2*pi)`，正负 kappa 均参与。")
    println(io,"\n## 单波数响应与完整局部响应")
    println(io,"- [roughness_spectrum.csv](roughness_spectrum.csv) 中的 `G_single` 是单位粗糙度 Fourier 分量的速度能量响应，`abs_Hhat_sq_G_single` 是 Gaussian 谱权重后的传递谱，不是某个固定 r 处的局部能量。")
    println(io,"- [reconstructed_response.csv](reconstructed_response.csv) 是先保留复数相位、再按 r 重构的完整响应；[local_energy.csv](local_energy.csv) 是重构后再计算的 `E_hf(r)`。")
    println(io,"- 物理径向积分能量：`$(main.Ephys)`；离散 Parseval 能量：`$(main.Eparse)`；相对差：`$(abs(main.Ephys-main.Eparse)/main.Eparse)`。")
    println(io,"- 以同一截断 Fourier 谱构造的壁面强迫一致性最大相对误差：`$(maximum(main.wall_error))`，比较的是重构 `q_hf(r,z=0)` 与 `-h_truncated(r) q_f,z(0)`。截断谱相对于解析 Gaussian 的最大几何误差为 `$(maximum(main.geometry_error))`。")
    println(io,"\n## 重新解释单 kappa 峰值")
    println(io,"此前 `|kappa|≈0.185` 的结果是单个单位 Fourier 粗糙度分量的传递响应峰值。完整 Gaussian 粗糙元的谱贡献由 `|Hhat(kappa)|^2 G_single(kappa)` 决定；因此峰值还要与 Gaussian 频谱宽度共同作用，不能直接把 `G_single` 的峰值当成局部粗糙元总能量峰值。")
    println(io,"负、正 kappa 响应在固定复数入射涡波下不必共轭；实值 Gaussian 的几何谱具有共轭对称性，但耦合响应还包含固定 q_f、目标算子和相位因子。")
    println(io,"\n## 重构收敛")
    println(io,"见 [reconstruction_convergence.csv](reconstruction_convergence.csv)。计算采用均匀 kappa 网格，`L_r=2*pi/Delta-kappa` 的周期径向区间，因此物理径向积分与离散 Parseval 比较使用同一周期区间。")
end
println("wrote Gaussian reconstruction results")

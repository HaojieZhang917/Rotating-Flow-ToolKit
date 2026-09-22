module CompressibleBEKBaseflowResiduals

using LinearAlgebra

export baseflow_primitive_residuals, residual_metrics

function residual_metrics(residual,term_scale,mask)
    selected=collect(mask)
    absolute=maximum(abs.(residual[selected]))
    scale=maximum(abs.(term_scale[selected]))
    term_relative=absolute/max(scale,eps(Float64))
    # A unit floor prevents a roundoff-size, one-term equation from being
    # reported as O(1) merely because residual and scale vanish together.
    relative=absolute/max(scale,1.0)
    rms=sqrt(sum(abs2,residual[selected])/length(selected))
    return (absolute_linf=absolute,relative_linf=relative,
        term_relative_linf=term_relative,rms=rms,scale_linf=scale)
end

"""Evaluate transformed and pre-DH primitive BEK base-flow residuals.

Inputs `U,V,W` use Lingwood's direct variables, scaled by the signed
differential rotation, on the transformed coordinate `y`. `H=T*W` is the
corresponding pre-DH normal-velocity profile; its velocity in the system-
rotation scale is `Ro*H`. `f,q` are the two thermal decomposition functions.
The supplied derivative matrices must act on the same `y` grid.
"""
function baseflow_primitive_residuals(y,D,D2,U,V,W,T,rho,f,q;
        Ro,Mr,Tw,gamma=1.4,Pr=0.72,max_y=12.0)
    n=length(y)
    all(length(v)==n for v in (U,V,W,T,rho,f,q)) ||
        throw(DimensionMismatch("all profiles must have the same length as y"))
    size(D)==(n,n) && size(D2)==(n,n) ||
        throw(DimensionMismatch("derivative matrices do not match the grid"))

    Co=2-Ro-Ro^2
    DU,DV,DW,DT=D*U,D*V,D*W,D*T
    D2U,D2V,D2T=D2*U,D2*V,D2*T
    Df,Dq,D2f,D2q=D*f,D*q,D2*f,D2*q
    shear=DU .^ 2 .+ DV .^ 2
    a=(gamma-1)*Mr^2/2
    radial_T=-2a.*f

    transformed_continuity=2 .* U .+ DW
    transformed_radial=Ro .* (U .^ 2 .+ W .* DU .- (V .^ 2 .- 1)) .-
        Co .* (V .- 1) .- D2U
    transformed_azimuthal=Ro .* (2 .* U .* V .+ W .* DV) .+ Co .* U .- D2V
    literature_f=D2f .- Pr*Ro .* W .* Df .- 2Pr*Ro .* U .* f .- 2Pr .* shear
    literature_q=D2q .- Pr*Ro .* W .* Dq
    current_f=D2f .+ Pr*Ro .* W .* Df .+ 2Pr*Ro .* U .* f .- 2Pr .* shear
    current_q=D2q .+ Pr*Ro .* W .* Dq
    literature_energy=D2T .- Pr*Ro .* W .* DT .- Pr*Ro .* U .* radial_T .+
        (gamma-1)*Pr*Mr^2 .* shear
    current_energy=D2T .+ Pr*Ro .* W .* DT .+ Pr*Ro .* U .* radial_T .+
        (gamma-1)*Pr*Mr^2 .* shear

    deta(v)=rho.*(D*v)
    H=T .* W
    mu=T
    k=T./Pr
    detaU,detaV,detaT=deta(U),deta(V),deta(T)
    primitive_continuity=2rho.*U.+deta(rho.*H)
    radial_inertia=rho .* (Ro .* (U .^ 2 .+ H .* detaU .- (V .^ 2 .- 1)) .-
        Co .* (V .- 1))
    radial_viscous=deta(mu .* detaU)
    primitive_radial=radial_inertia .- radial_viscous
    azimuthal_inertia=rho .* (Ro .* (2 .* U .* V .+ H .* detaV) .+ Co .* U)
    azimuthal_viscous=deta(mu .* detaV)
    primitive_azimuthal=azimuthal_inertia .- azimuthal_viscous
    primitive_state=rho .* T .- 1
    energy_convection=rho .* Ro .* (H .* detaT .+ U .* radial_T)
    energy_conduction=deta(k .* detaT)
    energy_dissipation=(gamma-1)*Mr^2 .* mu .* (detaU .^ 2 .+ detaV .^ 2)
    primitive_energy=energy_convection .- energy_conduction .- energy_dissipation

    mask=findall(i->i>1 && i<n && isfinite(y[i]) && y[i]<=max_y,eachindex(y))
    isempty(mask) && error("residual mask is empty")
    scales=(
        transformed_continuity=abs.(2 .* U) .+ abs.(DW),
        transformed_radial=abs.(Ro .* (U .^ 2 .+ W .* DU .- (V .^ 2 .- 1))) .+
            abs.(Co .* (V .- 1)) .+ abs.(D2U),
        transformed_azimuthal=abs.(Ro .* (2 .* U .* V .+ W .* DV)) .+
            abs.(Co .* U) .+ abs.(D2V),
        literature_f=abs.(D2f) .+ abs.(Pr*Ro .* W .* Df) .+
            abs.(2Pr*Ro .* U .* f) .+ abs.(2Pr .* shear),
        literature_q=abs.(D2q) .+ abs.(Pr*Ro .* W .* Dq),
        current_f=abs.(D2f) .+ abs.(Pr*Ro .* W .* Df) .+
            abs.(2Pr*Ro .* U .* f) .+ abs.(2Pr .* shear),
        current_q=abs.(D2q) .+ abs.(Pr*Ro .* W .* Dq),
        literature_energy=abs.(D2T) .+ abs.(Pr*Ro .* W .* DT) .+
            abs.(Pr*Ro .* U .* radial_T) .+ abs.((gamma-1)*Pr*Mr^2 .* shear),
        current_energy=abs.(D2T) .+ abs.(Pr*Ro .* W .* DT) .+
            abs.(Pr*Ro .* U .* radial_T) .+ abs.((gamma-1)*Pr*Mr^2 .* shear),
        primitive_continuity=abs.(2 .* rho .* U) .+ abs.(deta(rho .* H)),
        primitive_radial=abs.(radial_inertia) .+ abs.(radial_viscous),
        primitive_azimuthal=abs.(azimuthal_inertia) .+ abs.(azimuthal_viscous),
        primitive_state=abs.(rho .* T) .+ 1,
        primitive_energy=abs.(energy_convection) .+ abs.(energy_conduction) .+
            abs.(energy_dissipation),
    )
    residuals=(
        transformed_continuity=transformed_continuity,
        transformed_radial=transformed_radial,
        transformed_azimuthal=transformed_azimuthal,
        literature_f=literature_f,literature_q=literature_q,
        current_f=current_f,current_q=current_q,
        literature_energy=literature_energy,current_energy=current_energy,
        primitive_continuity=primitive_continuity,
        primitive_radial=primitive_radial,
        primitive_azimuthal=primitive_azimuthal,
        primitive_state=primitive_state,
        primitive_energy=primitive_energy,
    )
    metrics=NamedTuple{keys(residuals)}(map(name->
        residual_metrics(getfield(residuals,name),getfield(scales,name),mask),keys(residuals)))
    return (residuals=residuals,metrics=metrics,mask=mask,Co=Co,
        max_y=max_y,H=H,radial_T=radial_T)
end

end

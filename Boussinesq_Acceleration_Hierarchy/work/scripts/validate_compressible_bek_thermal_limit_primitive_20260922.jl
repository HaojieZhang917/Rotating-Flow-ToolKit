using BSplineKit
using Printf
using LinearAlgebra

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "work", "src", "CompressibleBEKLinearStability.jl"))
include(joinpath(ROOT, "work", "src", "CompressibleBEKBaseflowResiduals.jl"))
const CBEK = CompressibleBEKLinearStability
const BRES = CompressibleBEKBaseflowResiduals
const OUTDIR = joinpath(ROOT,"work","results",
    "compressible_bek_thermal_limit_gate_20260922")
mkpath(OUTDIR)

function thermal_on_grid(source,y)
    fs = BSplineKit.interpolate(source.y,source.f,BSplineOrder(4))
    qs = BSplineKit.interpolate(source.y,source.q,BSplineOrder(4))
    f = [fs(point) for point in y]
    q = [qs(point) for point in y]
    f[1]=0.0; f[end]=0.0; q[1]=1.0; q[end]=0.0
    return f,q
end

rows = NamedTuple[]
for (label,Ro,Mr,L,N) in (("von_Karman",-1.0,0.1,40.0,120),
                            ("matched_Ro_m0p1",-0.1,0.03,160.0,240))
    p = CBEK.PhysicalParameters(Ro=Ro,Mr=Mr,Tw=1.0,R=116.3,
        beta=0.137,omega=0.0,gamma=1.4,Pr=0.72)
    n = CBEK.NumericalParameters(N=N,domain=:finite,ymax=L)
    c = CBEK.prepare_case(p,n;allow_unverified_ro=true)
    source = CBEK.source_profiles(Ro,p.Pr;ymax=L)
    f,q = thermal_on_grid(source,c.y)
    check = BRES.baseflow_primitive_residuals(c.y,c.D,c.D2,c.F,c.G,
        c.H./c.T,c.T,c.rho,f,q;Ro=Ro,Mr=Mr,Tw=1.0,
        gamma=p.gamma,Pr=p.Pr,max_y=12.0)
    for name in (:primitive_continuity,:primitive_radial,:primitive_azimuthal,
                 :primitive_state,:primitive_energy,:literature_f)
        metric = getfield(check.metrics,name)
        push!(rows,(label=label,Ro=Ro,Mr=Mr,L=L,N=N,
            equation=String(name),absolute_linf=metric.absolute_linf,
            term_relative_linf=metric.term_relative_linf))
        @printf("%-16s %-22s abs=%.3e term_ratio=%.3e\n",
            label,string(name),metric.absolute_linf,metric.term_relative_linf)
    end
end
open(joinpath(OUTDIR,"primitive_residuals.csv"),"w") do io
    println(io,join(keys(first(rows)),','))
    foreach(r->println(io,join(values(r),',')),rows)
end

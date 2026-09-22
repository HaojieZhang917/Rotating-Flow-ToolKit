using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "work", "results", "learning_manual_audit_20260907")
const PR = 0.72

function rhs(g, u, r)
    p, h = u
    s = (2-r-r*r)/2
    q = s+r*g
    chi = (s+r)^2/q^2
    dp = 3r*p/q-PR*r*h
    dh = (3r*p/q+r*(chi-PR)*h)/(chi*q)
    [dp, dh]
end

function endpoint(r, n)
    u = [1.0, 0.0]
    dg = 1/n
    for k in 0:n-1
        g = k*dg
        k1 = rhs(g,u,r)
        k2 = rhs(g+dg/2,u+dg*k1/2,r)
        k3 = rhs(g+dg/2,u+dg*k2/2,r)
        k4 = rhs(g+dg,u+dg*k3,r)
        u += dg*(k1+2k2+2k3+k4)/6
    end
    u
end

function solve(n)
    lo, hi = -0.48, -0.46
    flo = endpoint(lo,n)[1]
    @assert flo*endpoint(hi,n)[1] < 0
    for _ in 1:48
        mid = (lo+hi)/2
        fm = endpoint(mid,n)[1]
        if fm*flo > 0
            lo, flo = mid, fm
        else
            hi = mid
        end
    end
    r = (lo+hi)/2
    p1,h1 = endpoint(r,n)
    p0 = -1/h1
    s = (2-r-r*r)/2
    omega = s+r
    sp = -0.5-r
    slope = -2*(omega/s)*((sp+1)*s-omega*sp)/s^2
    tw = 2-(omega/s)^2
    r,tw,p0,slope,p0*p1
end

mkpath(OUT)
open(joinpath(OUT,"outer_independent.csv"),"w") do io
    println(io,"g_steps,Ro_t,Tw_t,p0,dTcurve_dRo,p_at_g1")
    for n in (250,500,1000,2000)
        values = solve(n)
        @printf(io,"%d,%.14e,%.14e,%.14e,%.14e,%.14e\n",n,values...)
        @printf("n=%d Ro=%.12f Tw=%.12f p0=%.12f slope=%.12f p1=%.3e\n",n,values...)
    end
end
open(joinpath(OUT,"arithmetic_checks.txt"),"w") do io
    dr20,dr10 = -0.3121655358706,-0.3103311653013
    dt20,dt10 = 0.3254373350358,0.3202758274503
    println(io,"Two-point intercepts (arithmetic only; NOT converged asymptotic coefficients):")
    println(io,"Ro intercept = ",2dr10-dr20)
    println(io,"Tw intercept = ",2dt10-dt20)
    slope = solve(1000)[4]
    println(io,"Relative slope discrepancy epsilon=.02 = ",abs((dt20/dr20)/slope-1))
    println(io,"Relative slope discrepancy epsilon=.01 = ",abs((dt10/dr10)/slope-1))
    println(io,"Mapping eta_last/N^2 = ",8/(1.4*pi^2))
    println(io,"At outer endpoint gamma*(Tw-1) = ",solve(1000)[2]-1)
end

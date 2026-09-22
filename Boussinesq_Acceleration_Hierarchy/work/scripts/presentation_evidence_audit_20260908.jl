#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "src", "BEKIsothermal.jl"))
include(joinpath(@__DIR__, "..", "src", "BEKConsistent.jl"))
using .BEKIsothermal, .BEKConsistent
using LinearAlgebra, Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "work", "results", "presentation_evidence_audit_20260908")
const PR = 0.72
const GAMMA = 1.0

sfun(r) = (2-r-r^2)/2
omegafun(r) = sfun(r)+r
tcurve(r) = 1 + (1-(omegafun(r)/sfun(r))^2)/GAMMA

function tcurve_derivative(r)
    s = sfun(r); omega = omegafun(r); sp = -0.5-r; op = sp+1
    -2*(omega/s)*((op*s-omega*sp)/s^2)/GAMMA
end

function rhs_y(u,r)
    g,p,h = u
    s = sfun(r); omega = omegafun(r); q = s+r*g; chi = omega^2/q^2
    pp = 3r*p^2/q-PR*r*h*p
    hp = (pp+chi*r*h*p)/(chi*q)
    [p,pp,hp]
end

function rk4_y(p0,r,ymax,dy)
    n = ceil(Int,ymax/dy); step = ymax/n; u = [0.0,p0,0.0]
    for _ in 1:n
        k1=rhs_y(u,r); k2=rhs_y(u+step*k1/2,r)
        k3=rhs_y(u+step*k2/2,r); k4=rhs_y(u+step*k3,r)
        u += step*(k1+2k2+2k3+k4)/6
    end
    u
end

residual_y(z,ymax,dy) = begin
    u=rk4_y(z[1],z[2],ymax,dy); [u[1]-1,u[3]+1]
end

function jacobian_y(z,ymax,dy; relstep=2e-6)
    J=zeros(2,2)
    for j in 1:2
        d=relstep*max(abs(z[j]),0.1); zp=copy(z); zm=copy(z)
        zp[j]+=d; zm[j]-=d
        J[:,j]=(residual_y(zp,ymax,dy)-residual_y(zm,ymax,dy))/(2d)
    end
    J
end

function solve_y(initial,ymax,dy)
    z=copy(initial)
    for it in 1:20
        R=residual_y(z,ymax,dy); norm(R,Inf)<2e-11 && return z,R,it
        dz=-(jacobian_y(z,ymax,dy)\R); damping=1.0
        while damping>=2.0^-16
            trial=z+damping*dz
            if norm(residual_y(trial,ymax,dy),Inf)<norm(R,Inf)
                z=trial; break
            end
            damping/=2
        end
    end
    error("outer Y shooting did not converge")
end

function rhs_g(g,u,r)
    p,h=u; s=sfun(r); omega=omegafun(r); q=s+r*g; chi=omega^2/q^2
    [3r*p/q-PR*r*h,
     (3r*p/q+r*(chi-PR)*h)/(chi*q)]
end

function endpoint_g(r,n)
    u=[1.0,0.0]; dg=1/n
    for k in 0:n-1
        g=k*dg; k1=rhs_g(g,u,r); k2=rhs_g(g+dg/2,u+dg*k1/2,r)
        k3=rhs_g(g+dg/2,u+dg*k2/2,r); k4=rhs_g(g+dg,u+dg*k3,r)
        u += dg*(k1+2k2+2k3+k4)/6
    end
    u
end

function solve_g(n)
    lo,hi=-0.48,-0.46; flo=endpoint_g(lo,n)[1]
    @assert flo*endpoint_g(hi,n)[1]<0
    for _ in 1:52
        mid=(lo+hi)/2; fm=endpoint_g(mid,n)[1]
        if fm*flo>0; lo,flo=mid,fm; else; hi=mid; end
    end
    r=(lo+hi)/2; p1,h1=endpoint_g(r,n)
    (r=r,Tw=tcurve(r),p0=-1/h1,p_at_1=p1)
end

function derivative_checks()
    sol=solve_consistent_isothermal(-0.5;degree=40,a=2.0,b=0.6,c=0.5,
        gamma=GAMMA,Pr=PR,tolerance=1e-10)
    sol=solve_consistent_fixed_tw(1.2,sol;tolerance=1e-10)
    state=BEKConsistent.pack(sol.fields); m=length(state)
    vector=normalize(sin.(collect(1.0:m))); direction=normalize(cos.(collect(1.0:m)))
    K=BEKConsistent._jacobian_directional_derivative(state,sol.Ro,sol.Tw,
        sol.gamma,sol.Pr,sol.operators,vector)
    analytic_h=K*direction
    analytic_r=BEKConsistent._residual_ro_derivative(state,sol.Ro,sol.Tw,
        sol.gamma,sol.Pr,sol.operators)
    open(joinpath(OUT,"derivative_checks.csv"),"w") do io
        println(io,"step,hessian_relative_error,Ro_derivative_relative_error")
        for step in (1e-4,3e-5,1e-5,3e-6,1e-6)
            _,Jp=BEKConsistent._residual_jacobian(state+step*direction,sol.Ro,
                sol.Tw,sol.gamma,sol.Pr,sol.operators)
            _,Jm=BEKConsistent._residual_jacobian(state-step*direction,sol.Ro,
                sol.Tw,sol.gamma,sol.Pr,sol.operators)
            finite_h=((Jp-Jm)*vector)/(2step)
            rp,_=BEKConsistent._residual_jacobian(state,sol.Ro+step,sol.Tw,
                sol.gamma,sol.Pr,sol.operators)
            rm,_=BEKConsistent._residual_jacobian(state,sol.Ro-step,sol.Tw,
                sol.gamma,sol.Pr,sol.operators)
            finite_r=(rp-rm)/(2step)
            eh=norm(analytic_h-finite_h)/norm(finite_h)
            er=norm(analytic_r-finite_r)/norm(finite_r)
            @printf(io,"%.12e,%.12e,%.12e\n",step,eh,er)
        end
    end
end

function main()
    mkpath(OUT)
    configs=[(20.0,0.01),(30.0,0.01),(40.0,0.005),(60.0,0.005),(60.0,0.0025)]
    initial=[0.75,-0.48]
    open(joinpath(OUT,"outer_y_shooting.csv"),"w") do io
        println(io,"Ymax,dY,p0,Ro_t,Tw_t,residual,sigma_min,sigma_max,max_jacobian")
        for (ymax,dy) in configs
            z,R,_=solve_y(initial,ymax,dy); J=jacobian_y(z,ymax,dy); sv=svdvals(J)
            @printf(io,"%.8f,%.8f,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e,%.14e\n",
                ymax,dy,z[1],z[2],tcurve(z[2]),norm(R,Inf),sv[end],sv[1],maximum(abs,J))
            initial=z
        end
    end
    open(joinpath(OUT,"outer_g_coordinate.csv"),"w") do io
        println(io,"g_steps,Ro_t,Tw_t,p0,p_at_g1,dTcurve_dRo")
        for n in (250,500,1000,2000,4000)
            x=solve_g(n)
            @printf(io,"%d,%.14e,%.14e,%.14e,%.14e,%.14e\n",
                n,x.r,x.Tw,x.p0,x.p_at_1,tcurve_derivative(x.r))
        end
    end
    derivative_checks()
    y=solve_y([0.654,-0.470],60.0,0.0025)[1]; g=solve_g(4000)
    open(joinpath(OUT,"summary.txt"),"w") do io
        @printf(io,"Y-shooting: Ro_t=%.14e Tw_t=%.14e p0=%.14e\n",y[2],tcurve(y[2]),y[1])
        @printf(io,"g-coordinate: Ro_t=%.14e Tw_t=%.14e p0=%.14e\n",g.r,g.Tw,g.p0)
        @printf(io,"cross-method differences: dRo=%.3e dTw=%.3e dp0=%.3e\n",
            y[2]-g.r,tcurve(y[2])-g.Tw,y[1]-g.p0)
        @printf(io,"analytic endpoint slope dTw/dRo=%.14e\n",tcurve_derivative(g.r))
    end
end

main()

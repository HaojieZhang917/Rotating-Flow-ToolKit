#!/usr/bin/env julia

using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUT = joinpath(ROOT, "work", "results", "generic_outer_family_20260909")
const PR = 0.72
const GAMMA = 1.0

sfun(r) = (2-r-r^2)/2
omegafun(r) = sfun(r)+r
tcurve(r) = 1 + (1-(omegafun(r)/sfun(r))^2)/GAMMA

function rhs_g(g,u,r)
    p,h=u
    s=sfun(r); omega=omegafun(r); q=s+r*g; chi=omega^2/q^2
    [3r*p/q-PR*r*h,
     (3r*p/q+r*(chi-PR)*h)/(chi*q)]
end

function backward_g(r,n; profile=false)
    u=[0.0,-1.0]
    dg=-1/n
    data=profile ? Matrix{Float64}(undef,n+1,5) : nothing
    profile && (data[1,:].=(1.0,u[1],u[2],sfun(r)+r,1.0))
    for k in 1:n
        g=1+(k-1)*dg
        k1=rhs_g(g,u,r)
        k2=rhs_g(g+dg/2,u+dg*k1/2,r)
        k3=rhs_g(g+dg/2,u+dg*k2/2,r)
        k4=rhs_g(g+dg,u+dg*k3,r)
        u += dg*(k1+2k2+2k3+k4)/6
        if profile
            gn=1+k*dg; q=sfun(r)+r*gn; theta=1+(1-omegafun(r)^2/q^2)/GAMMA
            data[k+1,:].=(gn,u[1],u[2],q,theta)
        end
    end
    profile ? (p0=u[1],h0=u[2],data=data) : (p0=u[1],h0=u[2])
end

function rhs_y(u,r)
    g,p,h=u
    s=sfun(r); omega=omegafun(r); q=s+r*g; chi=omega^2/q^2
    pp=3r*p^2/q-PR*r*h*p
    hp=(pp+chi*r*h*p)/(chi*q)
    [p,pp,hp]
end

function endpoint_y(p0,h0,r,ymax,dy)
    n=ceil(Int,ymax/dy); step=ymax/n; u=[0.0,p0,h0]
    for _ in 1:n
        k1=rhs_y(u,r); k2=rhs_y(u+step*k1/2,r)
        k3=rhs_y(u+step*k2/2,r); k4=rhs_y(u+step*k3,r)
        u += step*(k1+2k2+2k3+k4)/6
    end
    u
end

function root_h0(n; lo=-0.49,hi=-0.45)
    flo=backward_g(lo,n).h0; fhi=backward_g(hi,n).h0
    flo*fhi<0 || error("h0 root not bracketed: $flo $fhi")
    for _ in 1:55
        mid=(lo+hi)/2; fm=backward_g(mid,n).h0
        if flo*fm>0; lo=mid; flo=fm; else; hi=mid; end
    end
    r=(lo+hi)/2; x=backward_g(r,n)
    (;r,p0=x.p0,h0=x.h0,Tw=tcurve(r))
end

function main()
    mkpath(OUT)
    rossby=sort(unique(vcat(collect(-.50:.01:-.30),
        [-.48,-.475,-.472,-.471,-.4705,-.4704,-.4703674163,
         -.4703,-.4702,-.4700,-.469,-.467,-.465])))
    open(joinpath(OUT,"family_convergence.csv"),"w") do io
        println(io,"g_steps,Ro,Tw0,p0,h0,Y60_g_error,Y60_h_error,Y60_p")
        for n in (1000,2000,4000,8000,16000), r in rossby
            x=backward_g(r,n)
            endpoint=endpoint_y(x.p0,x.h0,r,60.0,.005)
            @printf(io,"%d,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                n,r,tcurve(r),x.p0,x.h0,endpoint[1]-1,endpoint[3]+1,endpoint[2])
        end
    end
    open(joinpath(OUT,"h0_zero_convergence.csv"),"w") do io
        println(io,"g_steps,Ro_zero,Tw_zero,p0,h0")
        for n in (250,500,1000,2000,4000,8000,16000)
            x=root_h0(n)
            @printf(io,"%d,%.14e,%.14e,%.14e,%.14e\n",n,x.r,x.Tw,x.p0,x.h0)
        end
    end
    for r in (-.50,-.47,-.4703674163,-.46,-.40,-.30)
        x=backward_g(r,16000;profile=true)
        path=joinpath(OUT,"profile_Ro_$(replace(@sprintf("%.6f",abs(r)),"."=>"p")).csv")
        open(path,"w") do io
            println(io,"g,p,h,q,theta")
            for i in axes(x.data,1)
                @printf(io,"%.12e,%.12e,%.12e,%.12e,%.12e\n",x.data[i,:]...)
            end
        end
    end
    println("Output: $OUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

#!/usr/bin/env julia

using SparseArrays
using LinearAlgebra
using Printf

const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const OUTER=joinpath(ROOT,"work","results","generic_outer_family_20260909",
    "family_convergence.csv")
const OUT=joinpath(ROOT,"work","results","generic_inner_first_order_20260909")
const PR=0.72
const GAMMA=1.0

function csvrows(path)
    lines=filter(!isempty,strip.(readlines(path))); header=split(first(lines),',')
    [Dict(header .=> split(line,',')) for line in lines[2:end]]
end
num(row,key)=parse(Float64,row[key])
sfun(r)=(2-r-r^2)/2
omegafun(r)=sfun(r)+r

function solve_inner(r,p0,h0,L,deta)
    n=ceil(Int,L/deta); dz=L/n; eta=collect(0:n).*dz; nn=n+1
    iF(i)=i; iG(i)=nn+i; itau=2nn+1; neq=2nn+1
    rows=Int[]; cols=Int[]; vals=Float64[]; rhs=zeros(neq); row=0
    add(rr,cc,v)=(push!(rows,rr);push!(cols,cc);push!(vals,v))
    s=sfun(r); omega=omegafun(r); chi0=(omega/s)^2; Co=2-r-r^2
    a=2chi0*s; b=GAMMA*s^2/r; c=chi0*Co
    Btheta=2omega^2*r*p0/(GAMMA*s^3)
    for i in 2:nn-1
        row+=1
        add(row,iF(i-1),1/dz^2); add(row,iF(i),-2/dz^2); add(row,iF(i+1),1/dz^2)
        add(row,iG(i),-a); add(row,itau,b); rhs[row]=-b*Btheta*eta[i]
        row+=1
        add(row,iG(i-1),1/dz^2); add(row,iG(i),-2/dz^2); add(row,iG(i+1),1/dz^2)
        add(row,iF(i),c)
    end
    row+=1; add(row,iF(1),1.0)
    row+=1; add(row,iG(1),1.0)
    row+=1; add(row,iF(nn),1.0)
    row+=1
    add(row,iG(nn),3/(2dz)); add(row,iG(nn-1),-4/(2dz)); add(row,iG(nn-2),1/(2dz)); rhs[row]=p0
    row+=1
    for i in 1:nn
        add(row,iF(i),dz*((i==1 || i==nn) ? 0.5 : 1.0))
    end
    rhs[row]=-h0/2
    @assert row==neq
    x=sparse(rows,cols,vals,neq,neq)\rhs
    F=x[1:nn]; G=x[nn+1:2nn]; tau=x[itau]
    H=zeros(nn)
    for i in 2:nn
        H[i]=H[i-1]-dz*(F[i-1]+F[i])
    end
    D1=-GAMMA*s^2*tau
    residual=norm(sparse(rows,cols,vals,neq,neq)*x-rhs,Inf)
    (;eta,F,G,H,tau,D1,Btheta,residual,dz)
end

function main()
    mkpath(OUT)
    allrows=csvrows(OUTER)
    outer=filter(r->round(Int,num(r,"g_steps"))==16000,allrows)
    configs=[(12.0,.01),(16.0,.01),(20.0,.01),(24.0,.01),(24.0,.005)]
    open(joinpath(OUT,"convergence.csv"),"w") do io
        println(io,"Ro,L,deta,p0,h0,tau1,D1,Btheta,H1_end,F1_end,G1prime_target,residual")
        for row in outer
            r=num(row,"Ro"); p0=num(row,"p0"); h0=num(row,"h0")
            for (L,deta) in configs
                q=solve_inner(r,p0,h0,L,deta)
                @printf(io,"%.12e,%.8f,%.8f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                    r,L,q.dz,p0,h0,q.tau,q.D1,q.Btheta,q.H[end],q.F[end],p0,q.residual)
                if L==24.0 && deta==.005
                    path=joinpath(OUT,"profile_Ro_$(replace(@sprintf("%.6f",abs(r)),"."=>"p")).csv")
                    open(path,"w") do pio
                        println(pio,"eta,F1,H1,G1,Theta1")
                        for i in eachindex(q.eta)
                            @printf(pio,"%.12e,%.12e,%.12e,%.12e,%.12e\n",q.eta[i],q.F[i],q.H[i],q.G[i],q.tau+q.Btheta*q.eta[i])
                        end
                    end
                end
            end
        end
    end
    println("Output: $OUT")
end

abspath(PROGRAM_FILE)==abspath(@__FILE__) && main()

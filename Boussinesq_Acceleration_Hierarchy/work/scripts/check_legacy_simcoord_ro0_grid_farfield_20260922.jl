using LinearAlgebra
using NonlinearEigenproblems
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const SOURCE = normpath(joinpath(ROOT, "..", "compress", "BEK", "roughness", "CRD_STA.jl"))
const OUTDIR = joinpath(ROOT, "work", "results", "legacy_simcoord_ro0_convergence_20260922")
include(SOURCE)

# This reproduces CRD_BF.Cheb at (a,cap)=(6,20), while exposing its two
# distinct choices: rational-map scale and clipped profile-sampling position.
# The derivative matrices never incorporate the clipping.
function legacy_grid(N::Int, a::Real, cap::Real)
    theta = range(0, stop=pi, length=N+1)
    xi = reshape(-cos.(theta), N+1, 1)
    weights = [2; ones(N-1, 1); 2] .* (-1) .^ (0:N)
    X = repeat(xi, 1, N+1)
    dX = X - X'
    D = (weights * (1 ./ weights)') ./ (dX .+ I(N+1))
    D = D - diagm(vec(sum(D, dims=2)))
    b, c = 0.6, 0.5
    for i in 1:N+1
        s = xi[i]
        den = 1-b*s-(1-b)*(s^3+c*(1-s^2))
        D[i,:] .*= den^2 / (2a*(b+3*(1-b)*s^2-2c*(1-b)*s))
    end
    x = similar(xi)
    for i in 1:N+1
        s = xi[i]
        x[i] = min(a*(1+b*s+(1-b)*(s^3+c*(1-s^2))) /
            (1-b*s-(1-b)*(s^3+c*(1-s^2))), cap)
    end
    return D, D^2, x
end

function retained_indices(N)
    removed = (N+1,N+2,2N+2,2N+3,3N+3,3N+4,4N+4,4N+5,5N+5)
    setdiff(1:5*(N+1), removed)
end

function sampled_mode(eigvec, x, N; ymax=8.0, ns=161)
    full = zeros(ComplexF64, 5*(N+1))
    full[retained_indices(N)] = eigvec
    xx = vec(x)
    yy = range(0, ymax, length=ns)
    result = Vector{ComplexF64}(undef, 5*ns)
    for field in 0:4, (j,y) in enumerate(yy)
        hi = searchsortedfirst(xx, y)
        hi = clamp(hi, 2, N+1)
        lo = hi-1
        frac = (y-xx[lo])/(xx[hi]-xx[lo])
        result[field*ns+j] = (1-frac)*full[field*(N+1)+lo] +
                             frac*full[field*(N+1)+hi]
    end
    return result / norm(result)
end

function main()
    mkpath(OUTDIR)
    Ro, Co, Tw, Mr, R, beta, omega = 0.0, 2.0, 1.0, 0.1, 116.3, 0.137, 0.0
    gamma, sigma = 1.4, 0.72
    # One unchanged legacy BVP solution is reused in every stability run.
    u0,v0,w0,f,q,_,_,_ = baseflow_var(69,Ro,Co)
    Hsource,Tsource = T_ca(Mr,f,q,w0,gamma,Tw)
    cases = [("reference",69,6.0,20.0)]
    append!(cases, [("N",N,6.0,20.0) for N in (49,59,79,89)])
    append!(cases, [("map_scale",69,a,20.0) for a in (4.0,8.0,10.0)])
    append!(cases, [("sample_cap",69,6.0,cap) for cap in (30.0,40.0)])
    refprofile = nothing
    refalpha = 0.5403618717476425 + 0.00033737176837288875im
    rows = String[]
    for (axis,N,a,cap) in cases
        @printf("START %-11s N=%d a=%.1f cap=%.1f\n",axis,N,a,cap)
        flush(stdout)
        D,D2,x = legacy_grid(N,a,cap)
        if axis == "reference"
            Dold,D2old,xold = CRD_BF.Cheb(N)
            @printf("  legacy-grid reproduction: max|dD|=%.3e max|dD2|=%.3e max|dx|=%.3e\n",
                maximum(abs.(D-Dold)),maximum(abs.(D2-D2old)),maximum(abs.(x-xold)))
            @assert isapprox(D,Dold;rtol=1e-12) &&
                    isapprox(D2,D2old;rtol=1e-12) &&
                    isapprox(x,xold;rtol=1e-12)
        end
        F,G,H,T,rho,_ = interp(u0,v0,Hsource,Tsource,x,N,"sim")
        lam = -(2/3).*T
        kappa = T./sigma
        cof = Spatial_mode_BEK(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Mr/R,N,Ro,Co,D,D2)
        raw = assemble_mat(cof,D,D2,beta,omega)
        A = boudary_condition(raw...,N)
        vals,vecs = iar(PEP([A[1],A[2],A[3]]);
            σ=ComplexF64(0.5),neigs=4,maxit=500,tol=1e-12)
        profiles = [sampled_mode(vecs[:,j],x,N) for j in eachindex(vals)]
        if axis == "reference"
            pick = argmin(abs.(vals .- refalpha))
            refprofile = profiles[pick]
        else
            pick = argmax([abs(dot(refprofile,p)) for p in profiles])
        end
        alpha = vals[pick]
        overlap = abs(dot(refprofile,profiles[pick]))
        v = vecs[:,pick]
        backward = norm((A[1]+alpha*A[2]+alpha^2*A[3])*v) /
            ((opnorm(A[1])+abs(alpha)*opnorm(A[2])+abs2(alpha)*opnorm(A[3]))*norm(v))
        duplicates = count(==(cap),vec(x))
        row = join((axis,N,a,cap,real(alpha),imag(alpha),abs(alpha-refalpha),
                    overlap,backward,duplicates,length(vals)),',')
        push!(rows,row)
        open(joinpath(OUTDIR,"modes.csv"),"w") do io
            println(io,"axis,N,map_scale_a,sample_cap,alpha_real,alpha_imag,delta_from_notebook,overlap_with_reference,backward_error,cap_nodes,num_returned")
            foreach(r -> println(io,r), rows)
        end
        @printf("  alpha=%.12f%+.12fi overlap=%.8f backward=%.3e cap_nodes=%d\n",
            real(alpha),imag(alpha),overlap,backward,duplicates)
        flush(stdout)
    end
end

main()

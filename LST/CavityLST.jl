"""
CavityLST.jl
============

Pure-Julia port of the rotor–stator (two-disk, Batchelor-regime) local stability
model of Zhang & Liu, "Effects of surface mass flux on the stability and receptivity
of the rotor boundary layer in a rotor–stator system" (manuscript), following the
conventions of `TwoDisk/BaseFlow_cavity.jl` (mode 1) and `TwoDisk/Stability_Cavity.jl`:

  base flow   U = (r F(z), r G(z), H(z)) in units of Omega*delta, z = z*/h in [0, 1],
              Re_h = Omega h^2/nu = (h/delta)^2, rotor at z = 0, stator at z = 1;
                2F + H'/sqrt(Re_h) = 0,
                F''/Re_h - F' H/sqrt(Re_h) + G^2 - F^2 + Pi = 0,
                G''/Re_h - G' H/sqrt(Re_h) - 2 F G = 0,
              F(0) = 0, G(0) = 1, H(0) = -a_s;  F(1) = G(1) = H(1) = 0  (a_s > 0: suction).
              Solved here by Newton iteration on Chebyshev collocation with Pi as unknown
              (the original code uses the Pi-free differentiated form with scipy solve_bvp).
  operator    Spatial_mode_BEK1 evaluated with G - 1, i.e. in the ROTOR frame:
              omega_bar = R*omega = 0 is stationary on the rotor, the stator moves with V = -1.
              Dirichlet u = v = w = 0 on both disks, pressure gauge p = 0 at the stator:
              rows (1, N+1, N+2, 2N+2, 2N+3, 3N+3, 4N+4) removed (N = number of intervals).
  receptivity Gaussian wall displacement h = hr exp(-(r - r_f)^2/(2c^2)) e^{i(n theta - omega t)},
              transform Hx(alpha) = hr sqrt(pi/l_s) exp(-alpha^2/(4 l_s)), l_s = 1/(2c^2);
              wall velocity u_w = -F'(0) Hx, v_w = -G'(0) Hx, w_w = -i omega Hx (z-derivatives on [0,1]).

The residue is evaluated with the discrete transposed adjoint of the same constrained
pencil (Lambda = y^T b(alpha_n)/Q, Q = y^T (L1 + 2 alpha_n L2) x); with x normalised
max|v| = 1 this is the local coefficient |C_r| of the manuscript (Hermitian/continuous
adjoint there; the residue does not depend on the adjoint representation).
"""
module CavityLST

using LinearAlgebra

export cheb01, bary_interp, CavityBase, solve_cavity_baseflow, cavity_baseflow, sample_base,
       cavity_pencil, cavity_spectrum, wall_roughness_rhs, stator_eddy_viscosity, cavity_baseflow_turb, stator_sa_viscosity, residue_coefficient, full_state,
       interior_diff, wall_velocity_rhs, roughness_response, cavity_bilinear, cross_source, scattering_rhs

# ---------------------------------------------------------------- grid
"""Chebyshev–Gauss–Lobatto grid on [0,1] (z_1 = 0 rotor, z_{N+1} = 1 stator), D = d/dz, D2 = D^2."""
function cheb01(N::Int)
    x = -cos.(range(0, π; length=N+1))
    c = [2.0; ones(N-1); 2.0] .* (-1.0).^(0:N)
    X = repeat(x, 1, N+1); dX = X - X'
    Dx = (c*(1 ./ c)') ./ (dX + I)
    Dx -= Diagonal(vec(sum(Dx, dims=2)))
    z = (x .+ 1)/2; D = 2Dx
    return z, Matrix(D), Matrix(D*D)
end

"""Barycentric interpolation from a CGL grid on [0,1] with closed-form weights (overflow-free)."""
function bary_interp(z::AbstractVector, f::AbstractVector, zq::Real)
    N = length(z) - 1
    num = zero(eltype(f)); den = 0.0
    for j in 0:N
        d = zq - z[j+1]
        abs(d) < 1e-15 && return f[j+1]
        w = (-1.0)^j*(j == 0 || j == N ? 0.5 : 1.0)/d
        num += w*f[j+1]; den += w
    end
    return num/den
end

# ---------------------------------------------------------------- base flow
struct CavityBase
    Re_h::Float64
    a_s::Float64
    z::Vector{Float64}
    F::Vector{Float64}
    G::Vector{Float64}
    H::Vector{Float64}
    Pi::Float64
    D::Matrix{Float64}
    residual::Float64
    iterations::Int
end

"""
Newton–Chebyshev solution of the Batchelor similarity equations on N intervals.
`nut` (optional, length N+1): frozen eddy-viscosity ratio nu_t/nu at the nodes; the viscous terms
become [(1 + nu_t/nu) F']'/Re_h and [(1 + nu_t/nu) G']'/Re_h (used by `cavity_baseflow_turb`).
"""
function solve_cavity_baseflow(Re_h; a_s=0.0, N=300, guess=nothing, tol=1e-12, maxit=60, nut=nothing)
    z, D, D2 = cheb01(N); n = N + 1; s = sqrt(Re_h)
    mu = nut === nothing ? ones(n) : 1 .+ nut          # effective viscosity ratio
    dmu = nut === nothing ? zeros(n) : D*nut
    F, G, H, Pi = if guess === nothing
        Gc = 0.31
        G0 = Gc .+ (1 - Gc) .* exp.(-z .* s ./ 1.5) .- Gc .* exp.(-(1 .- z) .* s ./ 2.5)
        (zeros(n), G0, fill(-a_s, n) .* (1 .- z), -Gc^2)
    else
        (copy(guess.F), copy(guess.G), copy(guess.H), guess.Pi)
    end
    iF = 1:n; iG = n+1:2n; iH = 2n+1:3n; iP = 3n+1
    res = Inf; it = 0
    for k in 1:maxit
        it = k
        r = zeros(3n+1); J = zeros(3n+1, 3n+1)
        DF = D*F; DG = D*G; D2F = D2*F; D2G = D2*G; DH = D*H
        row = 0
        # radial momentum (interior)
        for i in 2:n-1
            row += 1
            r[row] = (mu[i]*D2F[i] + dmu[i]*DF[i])/Re_h - DF[i]*H[i]/s + G[i]^2 - F[i]^2 + Pi
            J[row, iF] .= (mu[i]*D2[i, :] .+ dmu[i]*D[i, :])/Re_h .- D[i, :]*H[i]/s; J[row, iF[i]] += -2F[i]
            J[row, iG[i]] = 2G[i]; J[row, iH[i]] = -DF[i]/s; J[row, iP] = 1
        end
        # azimuthal momentum (interior)
        for i in 2:n-1
            row += 1
            r[row] = (mu[i]*D2G[i] + dmu[i]*DG[i])/Re_h - DG[i]*H[i]/s - 2F[i]*G[i]
            J[row, iG] .= (mu[i]*D2[i, :] .+ dmu[i]*D[i, :])/Re_h .- D[i, :]*H[i]/s; J[row, iG[i]] += -2F[i]
            J[row, iF[i]] = -2G[i]; J[row, iH[i]] = -DG[i]/s
        end
        # continuity at all nodes but the stator; Pi absorbs the extra boundary condition on H
        for i in 1:n-1
            row += 1
            r[row] = DH[i] + 2s*F[i]
            J[row, iH] .= D[i, :]; J[row, iF[i]] = 2s
        end
        # boundary conditions
        for (idx, val) in ((iF[1], 0.0), (iF[n], 0.0), (iG[1], 1.0), (iG[n], 0.0), (iH[1], -a_s), (iH[n], 0.0))
            row += 1
            r[row] = (idx in iF ? F[idx - iF[1] + 1] : idx in iG ? G[idx - iG[1] + 1] : H[idx - iH[1] + 1]) - val
            J[row, idx] = 1
        end
        row == 3n+1 || error("row count $(row) ≠ $(3n+1)")
        res = norm(r, Inf)
        res < tol && break
        dx = J \ (-r)
        F .+= dx[iF]; G .+= dx[iG]; H .+= dx[iH]; Pi += dx[iP]
        # the residual floor is set by round-off in D2 (entries ~ N^4); stop on the Newton step instead
        norm(dx, Inf) < 1e-13 && (res = norm(r, Inf); break)
    end
    res < 1e-7 || error("base flow Newton did not converge (residual $(res))")
    return CavityBase(Float64(Re_h), Float64(a_s), z, F, G, H, Pi, D, res, it)
end

"""Base flow: direct Newton from the analytic guess; fine continuation in Re_h (factor 1.2) as fallback, then in a_s."""
function cavity_baseflow(Re_h; a_s=0.0, N=300)
    b = try
        solve_cavity_baseflow(Re_h; a_s=0.0, N=N)
    catch
        bb = solve_cavity_baseflow(100.0; N=N); Re = 100.0
        while Re < Re_h
            Re = min(1.2Re, Re_h); bb = solve_cavity_baseflow(Re; N=N, guess=bb)
        end
        bb
    end
    if a_s != 0
        for a in range(0, a_s, length=6)[2:end]
            b = solve_cavity_baseflow(Re_h; a_s=a, N=N, guess=b)
        end
    end
    return b
end

"""
Stator-side mixing-length eddy viscosity nu_t/nu on the nodes of base flow `b` at local radius R (in delta units).
Wall distance from the stator y = 1 - z (gap units); lab-frame shear magnitude S = sqrt(F'^2 + G'^2):
    nu_t/nu = (l/h)^2 R sqrt(Re_h) S * gamma(y),
    l/h     = lam delta_s tanh(kappa y Dvd/(lam delta_s))          (Michel et al. 1969 mixing length)
    Dvd     = 1 - exp(-y+/A+),  y+ = y sqrt(R sqrt(Re_h) S_w)        (van Driest damping, S_w = S at the stator)
    gamma   = 1/(1 + 5.5 (y/delta_s)^6)                             (Klebanoff intermittency; keeps rotor side laminar)
delta_s (outer length of the model): `edge = :first99` -> distance from the stator where G first reaches
0.99 K (K = G(z = 1/2)); `edge = :peak` -> distance of the Bödewadt overshoot maximum of G (the laminar
and weakly turbulent stator layer approaches the core with an overshoot, so :first99 marks only its inner part).
Default :peak (physical choice: the mixing-length outer scale should cover the whole layer including the overshoot). Note: lam = 0 gives nu_t = 0 (used to evaluate K, delta_s, S_w of any base flow).
"""
function stator_eddy_viscosity(b::CavityBase, R; lam=0.085, Aplus=26.0, kappa=0.41, edge=:peak)
    z = b.z; D = b.D; s = sqrt(b.Re_h)
    DF = D*b.F; DG = D*b.G; S = sqrt.(DF.^2 .+ DG.^2)
    K = bary_interp(z, b.G, 0.5)
    yq = range(0, 0.5, length=5001); Gq = [bary_interp(z, b.G, 1 - t) for t in yq]
    ds = if edge == :first99
        k = findfirst(g -> g >= 0.99K, Gq); k === nothing && error("stator-layer edge not found"); yq[k]
    elseif edge == :peak
        yq[argmax(Gq)]
    else
        error("edge must be :first99 or :peak")
    end
    y = 1 .- z; yplus = y .* sqrt(R*s*S[end])
    l = lam*ds .* tanh.(kappa .* y .* (1 .- exp.(-yplus ./ Aplus)) ./ (lam*ds))
    gam = 1 ./ (1 .+ 5.5 .* (y ./ ds).^6)
    return (nut=l.^2 .* (R*s) .* S .* gam, delta_s=ds, K=K, Sw=S[end])
end

"""
Spalart–Allmaras (SA-noft2) eddy viscosity for base flow `b` at local radius R, local similarity
(radial advection of the SA variable dropped; axial advection by w = Omega delta H kept).
chi = nu~/nu on [0,1] (gap units), wall distance to the stator d = 1 - z; in units nu^2/h^2:
  sqrt(Re_h) H chi' = c_b1 m(z) S~ chi - c_w1 f_w (chi/d)^2 + [((1+chi) chi')' + c_b2 chi'^2]/sigma,
  S~ = R sqrt(Re_h) sqrt(F'^2 + G'^2) + chi f_v2/(kappa d)^2   (boundary-layer shear, no solid-body vorticity),
chi(0) = chi(1) = 0. m(z) = (1 + tanh((z - z_m)/w_m))/2 switches production off in the rotor half
(laminar rotor layer, R below its transition radius). nu_t/nu = chi f_v1.
Pseudo-transient continuation (SER) with a finite-difference Jacobian; `chi0` is the initial guess.
"""
function stator_sa_viscosity(b::CavityBase, R; chi0=nothing, z_m=0.5, w_m=0.03, maxit=400, tol=1e-10, confine=true)
    cb1, sig, cb2, ka, cw2, cw3, cv1 = 0.1355, 2/3, 0.622, 0.41, 0.3, 2.0, 7.1
    cw1 = cb1/ka^2 + (1 + cb2)/sig
    z = b.z; D = b.D; n = length(z); s = sqrt(b.Re_h)
    DF = D*b.F; DG = D*b.G
    # boundary-layer form: shear magnitude |dV_h/dz| only. The axial vorticity 2 Omega G of the solid-body core
    # is excluded: standard SA produces nu~ in pure rotation (the reason for the SA-RC correction), which
    # would make the laminar core turbulent.
    Om = R*s .* sqrt.(DF.^2 .+ DG.^2)                                  # units nu/h^2
    m = (1 .+ tanh.((z .- z_m) ./ w_m)) ./ 2
    d = 1 .- z
    function resid(c)
        cp = max.(c, 0.0); dc = D*c
        r = zeros(n)
        fv1 = cp.^3 ./ (cp.^3 .+ cv1^3); fv2 = 1 .- cp ./ (1 .+ cp .* fv1)
        diff = D*((1 .+ cp) .* dc)
        for i in 2:n-1
            St = max(Om[i] + cp[i]*fv2[i]/(ka*d[i])^2, 0.3*Om[i])       # usual S~ clipping
            rr = St > 0 ? min(cp[i]/(St*(ka*d[i])^2), 10.0) : 10.0
            g = rr + cw2*(rr^6 - rr); fw = g*((1 + cw3^6)/(g^6 + cw3^6))^(1/6)
            r[i] = -s*b.H[i]*dc[i] + cb1*m[i]*St*cp[i] - cw1*fw*(cp[i]/d[i])^2 + (diff[i] + cb2*dc[i]^2)/sig
        end
        r[1] = c[1]; r[n] = c[n]
        return r
    end
    c = chi0 === nothing ? [y < 0.3 ? 20*y*(0.3 - y)/0.0225 : 0.0 for y in d] : copy(chi0)
    c[1] = 0; c[n] = 0
    # pseudo-transient continuation (switched evolution relaxation): (I/dt - J) dx = r, dt grows as r drops
    res = Inf; scale = max(1.0, maximum(abs.(Om))); dt = 1.0/scale; r0n = NaN
    for k in 1:maxit
        r0 = resid(c); rn = norm(r0, Inf); res = rn/scale
        res < tol && break
        isnan(r0n) && (r0n = rn)
        J = zeros(n, n)
        for j in 1:n
            h = 1e-7*max(1.0, abs(c[j])); cj = copy(c); cj[j] += h
            J[:, j] .= (resid(cj) .- r0) ./ h
        end
        M = -J; for i in 2:n-1; M[i, i] += 1/dt; end
        dx = M \ r0
        c .+= dx
        dt = min(1e12, (1.0/scale)*(r0n/max(rn, 1e-300))^1.5)
    end
    res < 1e-7 || error("SA Newton did not converge (residual $(res))")
    cp = max.(c, 0.0); nut = cp.^4 ./ (cp.^3 .+ cv1^3)
    if confine
        # the SA variable is hardly destroyed in the core and is carried towards the rotor by w < 0; without
        # confinement it turns the rotor layer turbulent (opposite K shift). Keep the rotor layer laminar:
        # Klebanoff intermittency on the stator-layer thickness (distance of the stator-side G maximum).
        yq = range(0, 0.5, length=5001); ds = yq[argmax([bary_interp(z, b.G, 1 - t) for t in yq])]
        nut = nut ./ (1 .+ 5.5 .* (d ./ ds).^6)
    end
    return (nut=nut, chi=c, residual=res)
end

"""
Base flow with a turbulent stator layer (laminar rotor layer and core), local similarity at radius R:
Picard iteration on the mixing-length eddy viscosity (`stator_eddy_viscosity`), Newton inside.
Starts from the laminar solution. Returns (base, nut, K, delta_s, delta_s_peak, delta_r, utau, delta_s_plus, iterations, change).
utau = friction velocity at the stator / (Omega r); delta_s_peak = distance of the stator-side G maximum;
delta_r = rotor-layer thickness = height of the first minimum of G above the rotor (Ekman undershoot).
"""
function cavity_baseflow_turb(Re_h, R; N=300, lam=0.085, Aplus=26.0, kappa=0.41, relax=0.5,
                              maxouter=400, tol=1e-10, guess=nothing, edge=:peak, model=:mixing_length,
                              chi0=nothing, confine=true)
    model in (:mixing_length, :sa) || error("model must be :mixing_length or :sa")
    b = guess === nothing ? cavity_baseflow(Re_h; N=N) : guess
    nut = zeros(N + 1); change = Inf; it = 0; chi = nothing
    if model == :sa        # seed SA from the converged mixing-length solution
        ml = cavity_baseflow_turb(Re_h, R; N=N, lam=lam, Aplus=Aplus, kappa=kappa, edge=edge)
        b = ml.base
        chi = chi0 !== nothing ? copy(chi0) :
              ml.nut .+ 10 .* (ml.nut .> 1e-3) .* (1 .- b.z) ./ max(ml.delta_s, 0.05) .* exp.(-(1 .- b.z) ./ max(ml.delta_s, 0.05))
    end
    for k in 1:maxouter
        it = k
        nn = if model == :mixing_length
            stator_eddy_viscosity(b, R; lam=lam, Aplus=Aplus, kappa=kappa, edge=edge).nut
        else
            sa = stator_sa_viscosity(b, R; chi0=chi, confine=confine); chi = sa.chi; sa.nut
        end
        new = (1 - relax) .* nut .+ relax .* nn
        change = norm(new - nut, Inf)/max(1.0, norm(new, Inf)); nut = new
        b = solve_cavity_baseflow(Re_h; N=N, guess=b, nut=nut)
        change < tol && break
    end
    change < 1e-8 || error("eddy-viscosity iteration did not converge (change $(change))")
    e = stator_eddy_viscosity(b, R; lam=lam, Aplus=Aplus, kappa=kappa, edge=edge)
    zq = range(0, 0.5, length=5001); Gq = [bary_interp(b.z, b.G, t) for t in zq]
    km = findfirst(k -> Gq[k] <= Gq[k-1] && Gq[k] <= Gq[k+1], 2:length(Gq)-1)
    dr = km === nothing ? NaN : zq[km+1]
    yq = range(0, 0.5, length=5001); dsp = yq[argmax([bary_interp(b.z, b.G, 1 - t) for t in yq])]
    s = sqrt(Re_h)
    return (base=b, nut=nut, K=e.K, delta_s=e.delta_s, delta_s_peak=dsp, delta_r=dr, utau=sqrt(e.Sw/(R*s)),
            delta_s_plus=e.delta_s*sqrt(R*s*e.Sw), iterations=it, change=change, chi=chi)
end

"""Sample F, G, H on another CGL grid on [0,1]."""
sample_base(b::CavityBase, zq) = ([bary_interp(b.z, b.F, x) for x in zq], [bary_interp(b.z, b.G, x) for x in zq],
                                  [bary_interp(b.z, b.H, x) for x in zq])

"""
d/dz on [0,1] for the degree N-2 polynomial through the N-1 interior CGL nodes (P_N-P_{N-2} pressure),
embedded in an (N+1)x(N+1) matrix with zero wall rows/columns. Closed-form barycentric weights
w_j = (-1)^j sin^2(j pi/N) (nodes are the zeros of U_{N-1}).
"""
function interior_diff(N::Int)
    x = -cos.((1:N-1)*π/N); w = [(-1.0)^j*sin(j*π/N)^2 for j in 1:N-1]; m = N - 1
    Dx = zeros(m, m)
    for i in 1:m, j in 1:m
        i == j && continue
        Dx[i, j] = (w[j]/w[i])/(x[i] - x[j])
    end
    for i in 1:m; Dx[i, i] = -sum(Dx[i, :]); end
    Dp = zeros(N+1, N+1); Dp[2:N, 2:N] .= 2Dx
    return Dp
end

# ---------------------------------------------------------------- operator (port of Spatial_mode_BEK1 + assemble_mat)
"""
Spatial pencil L0 + alpha L1 + alpha^2 L2 at local radius R, beta = n/R, omega (rotor frame).
`sigma` adds a complex frequency part (omega -> omega + i sigma) for Briggs–Bers tests.
`pressure = :PN` (default) reproduces the original code: pressure on all CGL nodes, continuity at all
nodes but the stator, gauge p(1) = 0. This P_N-P_N collocation carries a spurious zig-zag pressure mode
(p ~ T_N, whose derivative vanishes at interior nodes), harmless for eigenvalues but polluting forced
responses with weak core amplitudes. `pressure = :PN2` uses the P_N-P_{N-2} formulation: pressure only at
the interior nodes (degree N-2), continuity only at interior nodes, no gauge needed.
Returns raw 4(N+1) matrices, the kept (constrained) blocks, BC row indices, grid and base-flow samples.
"""
function cavity_pencil(b::CavityBase, R, beta, omega; N=199, sigma=0.0, pressure=:PN)
    z, D, D2 = cheb01(N); n = N + 1
    pressure in (:PN, :PN2) || error("pressure must be :PN or :PN2")
    F, Gl, H = sample_base(b, z); G = Gl .- 1.0          # rotor frame, as in the original code (G .- 1)
    s = sqrt(b.Re_h); Id = Matrix{Float64}(I, n, n); Z = zeros(n, n)
    dg(v) = Diagonal(v)
    DF = D*F; DG = D*G; DH = D*H
    Ta = [Id Z Z Z; Z Id Z Z; Z Z Id Z; Z Z Z Z]
    A = [dg(F) Z Z Id; Z dg(F) Z Z; Z Z dg(F) Z; Id Z Z Z]
    B = [dg(G/R) Z Z Z; Z dg(G/R) Z Id/R; Z Z dg(G/R) Z; Z Id/R Z Z]
    C = [dg(H/(R*s)) Z Z Z; Z dg(H/(R*s)) Z Z; Z Z dg(H/(R*s)) Id/s; Z Z Id/s Z]
    D1 = [dg(F/R) dg(-2(G .+ 1)/R) dg(DF/s) Z;
          dg(2(G .+ 1)/R) dg(F/R) dg(DG/s) Z;
          Z Z dg(DH/(R*s)) Z;
          Id/R Z Z Z]
    Vxx = [-Id/R Z Z Z; Z -Id/R Z Z; Z Z -Id/R Z; Z Z Z Z]
    Vyy = [-Id/R^3 Z Z Z; Z -Id/R^3 Z Z; Z Z -Id/R^3 Z; Z Z Z Z]
    Vzz = [-Id/(R*b.Re_h) Z Z Z; Z -Id/(R*b.Re_h) Z Z; Z Z -Id/(R*b.Re_h) Z; Z Z Z Z]
    K1 = kron(Matrix{Float64}(I, 4, 4), D); K2 = kron(Matrix{Float64}(I, 4, 4), D2)
    pressure == :PN2 && (K1[3n+1:4n, 3n+1:4n] .= interior_diff(N))   # dp/dz of the degree N-2 pressure
    om = omega + im*sigma
    L0 = ComplexF64.(D1) + im*R*beta*B - im*om*Ta - beta^2*R^2*Vyy + C*K1 + Vzz*K2
    L1 = ComplexF64.(im*A)
    L2 = ComplexF64.(-Vxx)
    bc = pressure == :PN ? [1, n, n+1, 2n, 2n+1, 3n, 4n] : [1, n, n+1, 2n, 2n+1, 3n, 3n+1, 4n]
    keep = setdiff(1:4n, bc)
    return (L0=L0[keep, keep], L1=L1[keep, keep], L2=L2[keep, keep], raw=(L0, L1, L2), keep=keep, bc=bc,
            z=z, D=D, F=F, G=Gl, H=H, n=n, N=N, R=R, beta=beta, omega=omega, Re_h=b.Re_h)
end

"""Full spectrum of the constrained quadratic pencil (companion linearisation), sorted by distance to `shift`."""
function cavity_spectrum(op; shift=0.3+0im, maxabs=10.0)
    m = size(op.L0, 1); Id = Matrix{ComplexF64}(I, m, m); Z = zeros(ComplexF64, m, m)
    E = eigen([-op.L1 -op.L0; Id Z], [op.L2 Z; Z Id])
    ids = [i for i in eachindex(E.values) if isfinite(E.values[i]) && abs(E.values[i]) < maxabs]
    sort!(ids, by=i -> abs(E.values[i] - shift))
    [(alpha=E.values[i], vector=E.vectors[m+1:2m, i]) for i in ids]
end

"""Embed a kept-row vector into the full 4(N+1) state (zeros on the constrained rows)."""
full_state(op, x) = (f = zeros(ComplexF64, 4op.n); f[op.keep] .= x; f)

"""
Right-hand side on the kept rows produced by the Gaussian wall displacement at complex alpha:
b = -P(alpha)[keep, bc] q_bc, q_bc = (u_w, v_w, w_w) at the rotor.
"""
function wall_roughness_rhs(op, alpha; c=1.0, hr=1/op.Re_h)
    ls = 1/(2c^2); Hx = hr*sqrt(π/ls)*exp(-alpha^2/(4ls))
    n = op.n; DF0 = (op.D*op.F)[1]; DG0 = (op.D*op.G)[1]
    q = zeros(ComplexF64, 4n); q[1] = -DF0*Hx; q[n+1] = -DG0*Hx; q[2n+1] = -im*op.omega*Hx
    P = op.raw[1] + alpha*op.raw[2] + alpha^2*op.raw[3]
    return -P[op.keep, op.bc]*q[op.bc]
end

"""Residue coefficient Lambda = y^T b / Q for direct (x) and left (y) eigenvectors at alpha_n."""
residue_coefficient(op, alpha, x, y, b) = (transpose(y)*b)/(transpose(y)*(op.L1 + 2alpha*op.L2)*x)

"""Kept-row right-hand side of a prescribed rotor-wall velocity (u_w, v_w, w_w) at alpha: -P[keep,bc] q_bc."""
function wall_velocity_rhs(op, alpha, uw, vw, ww=0.0)
    n = op.n; q = zeros(ComplexF64, 4n); q[1] = uw; q[n+1] = vw; q[2n+1] = ww
    P = op.raw[1] + alpha*op.raw[2] + alpha^2*op.raw[3]
    return -P[op.keep, op.bc]*q[op.bc]
end

"""
Local bilinear convection B(a, b) = (a.grad) b + curvature, for two local Fourier waves on the
same grid (full 4n states; perturbation velocities in local r*Omega units, as in the operator):
    B_r = a_u i alpha_b b_u + a_v i beta_b b_u + a_w D b_u/sqrt(Re_h) - a_v b_v/R
    B_t = a_u i alpha_b b_v + a_v i beta_b b_v + a_w D b_v/sqrt(Re_h) + a_v b_u/R
    B_z = a_u i alpha_b b_w + a_v i beta_b b_w + a_w D b_w/sqrt(Re_h)
(alpha_b, beta_b) are the wavenumbers of b. Coriolis is linear and does not appear.
With a = base flow (F, G-1, H/R) this reproduces the convective part of the operator exactly
(validate_cavity_scattering.jl, check V1).
"""
function cavity_bilinear(a, b, alpha_b, beta_b, op)
    n = op.n; s = sqrt(op.Re_h); R = op.R; D = op.D
    ua = a[1:n]; va = a[n+1:2n]; wa = a[2n+1:3n]
    ub = b[1:n]; vb = b[n+1:2n]; wb = b[2n+1:3n]
    adv = (im*alpha_b).*ua .+ (im*beta_b).*va
    Bu = adv.*ub .+ wa.*(D*ub)./s .- va.*vb./R
    Bv = adv.*vb .+ wa.*(D*vb)./s .+ va.*ub./R
    Bw = adv.*wb .+ wa.*(D*wb)./s
    return Bu, Bv, Bw
end

"""
First-order response to a unit rotor-steady roughness Fourier amplitude at complex `kappa` on `op`
(omega = 0): u(0) = -F'(0), v(0) = -G'(0), w(0) = 0, no-slip at the stator. Returns the FULL 4n field
including the wall values (they are part of the field and enter D q_h in the cross source).
"""
function roughness_response(op, kappa; h=1.0)
    n = op.n; DF0 = (op.D*op.F)[1]; DG0 = (op.D*op.G)[1]
    P = op.raw[1] + kappa*op.raw[2] + kappa^2*op.raw[3]
    q = zeros(ComplexF64, 4n)
    q[op.keep] .= P[op.keep, op.keep] \ wall_velocity_rhs(op, kappa, -h*DF0, -h*DG0, 0.0)
    q[1] = -h*DF0; q[n+1] = -h*DG0
    return q
end

"""Full 4n cross source -(B(q_h,q_f) + B(q_f,q_h)) of the (kappa, beta_h) and (alpha_f, beta_f) waves."""
function cross_source(op, qh, kappa, beta_h, qf, alpha_f, beta_f)
    N1 = cavity_bilinear(qh, qf, alpha_f, beta_f, op)
    N2 = cavity_bilinear(qf, qh, kappa, beta_h, op)
    return vcat(-(N1[1] .+ N2[1]), -(N1[2] .+ N2[2]), -(N1[3] .+ N2[3]), zeros(ComplexF64, op.n))
end

"""
Kept-row right-hand side of the second-order scattered field at (alpha_t = kappa + alpha_f,
beta_t = beta_h + beta_f, omega_t = omega_f) on the target operator `op`:
volume cross source plus the Taylor-expanded no-slip condition on the rough rotor,
q_hf(0) = -h dq_f/dz(0) (z in gap units, as in `wall_roughness_rhs`); q_hf(1) = 0.
`qh` is the response to a unit roughness Fourier amplitude; `h` scales the result.
Returns (total, volume, wall) so the two contributions can be separated.
"""
function scattering_rhs(op, alpha_t, qh, kappa, beta_h, qf, alpha_f, beta_f; h=1.0)
    n = op.n; D = op.D
    vol = h .* cross_source(op, qh, kappa, beta_h, qf, alpha_f, beta_f)[op.keep]
    du = (D*qf[1:n])[1]; dv = (D*qf[n+1:2n])[1]; dw = (D*qf[2n+1:3n])[1]
    wall = wall_velocity_rhs(op, alpha_t, -h*du, -h*dv, -h*dw)
    return (total=vol .+ wall, volume=vol, wall=wall)
end

end # module

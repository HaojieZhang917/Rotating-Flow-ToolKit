"""
local_mode_tracking.jl
======================

Local (shift-invert type) solver for one eigenvalue of the spatial quadratic
pencil L(alpha) = L0 + alpha*L1 + alpha^2*L2 from `VonKarmanLST.jl`.

Ruhe-type nonlinear inverse iteration: with L'(alpha) = L1 + 2*alpha*L2 and
a fixed normalisation vector e,

    u        = L(alpha_k) \\ (L'(alpha_k) x_k)
    alpha_k+1 = alpha_k - (e'x_k)/(e'u)
    x_k+1    = u/(e'u)

Each step costs one n×n LU of L(alpha), instead of a dense QZ of the 2n×2n
linearisation. Convergence is local (quadratic for simple eigenvalues), so it
is meant for refinement and continuation from a known eigenpair, not for
discovering new modes. Callers must keep the full-spectrum N-convergence
check for the starting mode.

Also provides:
- `continue_mode`: follow a mode along a parameter path with small steps;
- `bb_branch_local`: Briggs–Bers branch test, omega -> omega + i*sigma, using
  the local solver at each sigma (L0 gets +sigma*Ta on the kept rows).
"""

using LinearAlgebra

"""Refine an eigenpair of L0 + a*L1 + a^2*L2 near `alpha0`.

Returns (alpha, vector, residual, iterations, converged).
"""
function refine_mode(L0, L1, L2, alpha0::Number; x0=nothing, tol=1e-12, maxit=30)
    n = size(L0, 1)
    a = ComplexF64(alpha0)
    x = x0 === nothing ? ones(ComplexF64, n) : ComplexF64.(x0)
    if x0 === nothing
        # one inverse step from a generic vector to get into the eigenspace
        x = lu(L0 + a*L1 + a^2*L2) \ x
    end
    e = conj.(x) ./ max(norm(x)^2, eps())
    x ./= dot(conj.(e), x) == 0 ? 1 : transpose(e)*x
    conv = false; it = 0
    for k in 1:maxit
        it = k
        L = L0 + a*L1 + a^2*L2
        F = lu(L; check=false)
        issuccess(F) || (conv = true; break)          # exact singularity: a is an eigenvalue
        u = F \ ((L1 + 2a*L2)*x)
        s = transpose(e)*u
        da = (transpose(e)*x)/s
        a -= da
        x = u ./ s
        if abs(da) <= tol*max(abs(a), 1.0)
            conv = true; break
        end
    end
    r = norm((L0 + a*L1 + a^2*L2)*x)/max(norm(L0)*norm(x), eps())
    return (alpha=a, vector=x, residual=r, iterations=it, converged=conv)
end

"""
Follow a mode along `params` (any iterable). `build(p)` must return a named
tuple with fields L0, L1, L2 for parameter value p. Steps where the local
solve fails to converge or jumps by more than `max_jump` stop the path.
"""
function continue_mode(build, params, alpha0, x0=nothing; max_jump=0.05, kwargs...)
    out = NamedTuple[]
    a, x = ComplexF64(alpha0), x0
    for p in params
        op = build(p)
        m = refine_mode(op.L0, op.L1, op.L2, a; x0=x, kwargs...)
        (m.converged && abs(m.alpha - a) <= max_jump) || break
        push!(out, (param=p, alpha=m.alpha, vector=m.vector, residual=m.residual))
        a, x = m.alpha, m.vector
    end
    return out
end

"""
Briggs–Bers test with the local solver. `build(sigma)` returns the operator at
complex frequency omega + i*sigma. Returns (branch, Im alpha at last sigma, last sigma).
"""
function bb_branch_local(build, alpha0, x0, sigmas; max_jump=0.05)
    path = continue_mode(build, sigmas, alpha0, x0; max_jump=max_jump)
    isempty(path) && return ("undetermined", imag(alpha0), 0.0)
    a = path[end].alpha; s = path[end].param
    dai = imag(a) - imag(alpha0)
    br = if s == sigmas[end]
        imag(a) > 0 ? "plus" : "minus"
    elseif imag(a) > 0 && dai > 0
        "plus"
    elseif imag(a) < 0 && dai < 0
        "minus"
    else
        "undetermined"
    end
    return (br, imag(a), s)
end

using Test
using LinearAlgebra
using Printf

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(ROOT, "work", "src", "CompressibleBEKLinearStability.jl"))
const CBEK = CompressibleBEKLinearStability
const OUTDIR = joinpath(ROOT, "work", "results", "compressible_bek_thermal_limit_gate_20260922")
mkpath(OUTDIR)

function parameters(;Ro,Mr,Tw=1.0,R=116.3,beta=0.137)
    CBEK.PhysicalParameters(Ro=Ro,Mr=Mr,Tw=Tw,R=R,beta=beta,
        omega=0.0,gamma=1.4,Pr=0.72)
end

function numerical(;N=49,target=0.53+0im,thermal_model=:halfline,ymax=40.0)
    CBEK.NumericalParameters(N=N,domain=:rational,ymax=ymax,
        target=ComplexF64(target),neigs=3,tol=1e-10,maxit=600,
        balance=true,balance_iterations=2,regularization=0.0,
        thermal_model=thermal_model)
end

function max_structure_error(gates)
    maximum(abs, (gates.continuity_curvature,
        gates.continuity_base_transport,gates.radial_coriolis,
        gates.azimuthal_coriolis,gates.energy_wallnormal_scale))
end

function mode_row(label,p,n)
    case = CBEK.prepare_case(p,n;allow_unverified_ro=true)
    gates = CBEK.operator_structure_gates(case)
    derivative = CBEK.directional_qep_test(case)
    result = CBEK.solve_mode(case)
    alpha = first(result.alpha)
    source = CBEK.source_profiles(p.Ro,p.Pr;ymax=n.ymax)
    row = (label=label,Ro=p.Ro,Mr=p.Mr,R=p.R,beta=p.beta,N=n.N,
        thermal_model=String(n.thermal_model),
        evidence_scope=String(case.evidence_scope),
        farfield_exponent=p.Pr*p.Ro*last(source.w),
        max_temperature_change=maximum(abs.(case.T .- 1)),
        alpha_real=real(alpha),alpha_imag=imag(alpha),
        qep_backward_error=first(result.backward_errors),
        operator_structure_error=max_structure_error(gates),
        qep_derivative_error=derivative.relative_error)
    @printf("%-19s Ro=%6.3f Mr=%.4f alpha=%.9f%+.9fi Tdev=%.3e scope=%s\n",
        label,p.Ro,p.Mr,real(alpha),imag(alpha),
        row.max_temperature_change,row.evidence_scope)
    @test row.operator_structure_error < 1e-8
    @test row.qep_derivative_error < 1e-5
    @test row.qep_backward_error < 1e-8
    return row
end

rows = NamedTuple[]
@testset "thermal half-line guards and diagnostic scope" begin
    for p in (parameters(Ro=0.0,Mr=0.1),
              parameters(Ro=0.0,Mr=0.0,Tw=1.1))
        @test_throws ArgumentError CBEK.prepare_case(p,numerical();
            allow_unverified_ro=true)
    end
    p0 = parameters(Ro=0.0,Mr=0.0)
    c0 = CBEK.prepare_case(p0,numerical();allow_unverified_ro=true)
    @test maximum(abs.(c0.T .- 1)) < 1e-12
    @test_throws ArgumentError CBEK.assemble_qep(c0)
    @test CBEK.at_fixed_system_mach(p0,0.0,0.3).Mr == 0.0
    @test isapprox(CBEK.at_fixed_system_mach(p0,-0.1,0.3).Mr,0.03)
    @test_throws ArgumentError CBEK.prepare_case(
        parameters(Ro=1.0,Mr=0.1),numerical();allow_unverified_ro=true)
    @test_throws ArgumentError CBEK.prepare_case(
        parameters(Ro=-0.1,Mr=0.03),numerical();allow_unverified_ro=true)
end

@testset "operator and matched-limit numerical gates" begin
    push!(rows,mode_row("von_Karman",parameters(Ro=-1.0,Mr=0.1),
        numerical(target=0.38+0im)))
    reference = parameters(Ro=0.0,Mr=0.0)
    for (Ro,L) in ((-0.2,80.0),(-0.1,160.0),(-0.05,320.0))
        p = CBEK.at_fixed_system_mach(reference,Ro,0.3)
        push!(rows,mode_row("matched_Msystem",p,numerical(ymax=L)))
    end
    push!(rows,mode_row("matched_Msystem_N59",
        CBEK.at_fixed_system_mach(reference,-0.05,0.3),
        numerical(N=59,ymax=320.0)))
    push!(rows,mode_row("Ro1_finite_diagnostic",
        parameters(Ro=1.0,Mr=0.1),
        numerical(target=0.49+0im,thermal_model=:finite_diagnostic)))
    @test rows[end].evidence_scope == "thermal_outflow_unverified"
    push!(rows,mode_row("Ro0_finite_diagnostic",
        parameters(Ro=0.0,Mr=0.1),
        numerical(N=69,target=0.54+0im,
            thermal_model=:finite_diagnostic)))
    @test rows[end].evidence_scope == "ro0_thermal_degenerate"
    @test abs((rows[end].alpha_real+im*rows[end].alpha_imag) -
              (0.530157206567-0.000748070374im)) < 1e-4
end

open(joinpath(OUTDIR,"cases.csv"),"w") do io
    println(io,join(keys(first(rows)),','))
    foreach(r->println(io,join(values(r),',')),rows)
end
println("Wrote ", OUTDIR)

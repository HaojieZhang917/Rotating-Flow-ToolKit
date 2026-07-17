module LegacyNeutralCurveEntry
include(joinpath(@__DIR__, "..", "NeutralCurveRunner.jl"))
using .NeutralCurveRunner
const config = CurveConfig(Tw=1.0)
end

module LegacySutherlandEntry
include(joinpath(@__DIR__, "..", "SutherlandMarching.jl"))
using .SutherlandMarching
const params = Params(verbose=false)
end

@testset "Legacy entry points" begin
    @test LegacyNeutralCurveEntry.config.model == :lopez
    @test LegacySutherlandEntry.params.Tw == 1.5
end

using Test
using RotatingDiskFlow

@testset "RotatingDiskFlow" begin
    include("test_api.jl")
    include("test_numerics.jl")
    include("test_lopez_operator.jl")
    include("test_compatibility.jl")

    if lowercase(get(ENV, "RUN_PHYSICS_REGRESSION", "false")) in
       ("1", "true", "yes")
        include("test_malik_benchmarks.jl")
    else
        @info "Skipping Malik benchmarks; set RUN_PHYSICS_REGRESSION=true to run them"
    end
end

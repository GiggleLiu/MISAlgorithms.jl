using MISAlgorithms
using Test

@testset "Core" begin
    include("Core.jl")
    include("EliminateGraph.jl")
    include("mis1.jl")
end

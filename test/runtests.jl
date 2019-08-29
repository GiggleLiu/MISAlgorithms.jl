using MISAlgorithms
using Test

@testset "Core" begin
    include("Core.jl")
end

@testset "EliminateGraph" begin
    include("EliminateGraph.jl")
end

@testset "mis1" begin
    include("mis1.jl")
end

@testset "graphlib" begin
    include("graphlib.jl")
end

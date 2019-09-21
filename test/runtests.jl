using MISAlgorithms
using Test

@testset "Core" begin
    include("Core.jl")
end

@testset "EliminateGraph" begin
    include("EliminateGraph.jl")
end

@testset "graphlib" begin
    include("graphlib.jl")
end

@testset "mis1" begin
    include("mis1.jl")
end

@testset "mis2" begin
    include("mis2.jl")
end

@testset "binary sparse tensor" begin
    include("SparseTensors/BinarySparseTensor.jl")
end

using MISAlgorithms
using Test

@testset "mis2" begin
    @test MISAlgorithms.minx(x->x^2, [2,-1, 3]) == -1

    graph = EliminateGraph([0 1 0 0 0;
                            1 0 1 1 0;
                            0 1 0 1 0;
                            0 1 1 0 1;
                            0 0 0 1 0])

    @test mis2(graph) == 3
    graph = disconnected_cliques_eg(6, 6)
    @test mis2(graph) == 2
end

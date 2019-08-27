using MISAlgorithms

@testset "mis1" begin
    @test MISAlgorithms.minx(x->x^2, [2,-1, 3]) == -1

    graph = EliminateGraph([0 1 0 0 0;
                            1 0 1 1 0;
                            0 1 0 1 0;
                            0 1 1 0 1;
                            0 0 0 1 0])

    @test mis1(graph) == 3
    #@test mis2(graph) == 3
end

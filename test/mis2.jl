using MISAlgorithms
using Test

@testset "mis2" begin
    @test MISAlgorithms.minx(x->x^2, [2,-1, 3]) == -1

    g1 = EliminateGraph([0 1 0 0 0;
                            1 0 1 1 0;
                            0 1 0 1 0;
                            0 1 1 0 1;
                            0 0 0 1 0])

    @test mis2(g1) == 3
    g2 = disconnected_cliques_eg(6, 6)
    @test mis2(g2) == 2
    g3 = K_eg(7, 2)
    @test mis2(g3) == 7
end

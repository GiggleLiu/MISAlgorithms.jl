using Test
using MISAlgorithms

@testset "constructors and properties" begin
    tbl = Bool[false true false false; true false true true; false true false true; false true true false]
    eg = EliminateGraph(tbl)
    @test check_validity(eg)
    @test eg == EliminateGraph(4,[1=>2, 2=>3,2=>4,3=>4])
    @test collect(Vertices(eg)) == collect(1:4)
end

@testset "neighbors" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @show eg
    nbs = neighbors(eg, 1)
    @test nbs == [2]
    @test neighbors2(eg, 1, nbs) == [3,4]
    @test neighborcover(eg, 1) == [1,2]
    @test nv(eg) == 4
end

@testset "degree" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @test degrees(eg) == [1, 3, 2, 2]

    vs = Vertices(eg)
    @test [vs...] == vertices(eg)
    @test vs[3] == vertices(eg)[3]
    vmin, dmin = mindegree_vertex(eg)
    @test degree(eg, vmin) == minimum(degrees(eg)) == dmin
    vmax, dmax = maxdegree_vertex(eg)
    @test degree(eg, vmax) == maximum(degrees(eg)) == dmax
    vmin, vmax, dmin, dmax = minmaxdegree_vertex(eg)
    @test degree(eg,vmin) == minimum(degrees(eg)) == dmin
    @test degree(eg,vmax) == maximum(degrees(eg)) == dmax
end

@testset "vertices and neighbors - eliminated" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    eliminate!(eg, 2)
    eliminate!(eg, 3)
    @test nv(eg) == 2
    recover!(eg)
    @test nv(eg) == 3
    eg2 = eg \ 3
    @test nv(eg) == 3
    @test nv(eg2) == 2
    @test vertices(eg2) == [1,4]
    @test neighbors(eg2, 1) == []
    @test degrees(eg) == [1, 0, 1]
    @test degrees(eg2) == [0, 0]
    eg3 = eg2 \ (1,4)
    @test nv(eg3) == 0
end

@testset "eliminate and recover" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @test check_validity(eg)
    res = eliminate(eg, 3) do eg
        @test check_validity(eg)
        nv(eg)
    end # do
    @test res == 3

    eg4 = eg \ NeighborCover(1)
    @test vertices(eg4) == [3,4]

    res = eliminate(eg, (3,4)) do eg
        @test check_validity(eg)
        eliminate(eg, (1,2)) do eg
            @test check_validity(eg)
            nv(eg)
        end
    end
    @test res == 0
    @test nv(eg) == 4
    @test eg.level == 0
end

@testset "mirrors" begin
    # isclique
    g = disconnected_cliques_eg(3,4)
    @test isclique(g, [1,2,3])
    @test !isclique(g, [1,2,4])
    @test !isclique(g)
    g = K_eg(3)
    @test isclique(g)

    nbs = neighbors(petersen_graph, 1)
    @test mirrors(petersen_graph, 1, nbs) == []
    @test mirrorcover(petersen_graph, 1) == [1]

    @test Set(mirrorcover(K_eg(3,3), 1)) == Set([1,2,3])
    @test Set(mirrorcover(K_eg(3,3), 4)) == Set([4,5,6])
end

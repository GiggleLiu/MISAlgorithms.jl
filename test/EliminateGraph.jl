using Test
using MISAlgorithms

@testset "constructors" begin
    eg = rand_eg(10, 0.5)
    @test eg isa EliminateGraph
    @test nv(eg) == nv0(eg) == 10
end

@testset "vertices" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @test neighbors(eg, 1) == [2]
    @test neighborcover(eg, 1) == [1,2]
    @test nv(eg) == 4
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
    res = eliminate(eg, 3) do eg
        nv(eg)
    end # do
    @test res == 3

    eg4 = eg \ NeighborCover(1)
    @test vertices(eg4) == [3,4]

    res = eliminate(eg, (3,4)) do eg
        eliminate(eg, (1,2)) do eg
            nv(eg)
        end
    end
    @test res == 0
    @test nv(eg) == 4
    @test eg.level == 0
end

@testset "degree" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')

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

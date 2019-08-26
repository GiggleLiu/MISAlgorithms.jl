include("utils.jl")
include("EliminateGraph.jl")

function mis1(eg::EliminateGraph)
    N = nv(eg)
    if N == 0
        return 0
    else
        v= mindegree_vertex(eg)
        return 1 + neighborcover_mapreduce(y->eliminate(mis1, eg, NeighborCover(y)), max, eg, v)
    end
end

function mindegree_vertex(eg::EliminateGraph)
    N0 = nv0(eg)
    N = nv(eg)
    dmin = 999999
    vmin = 0
    for j = N0-N+1:N0
        @inbounds vj = eg.vertices[j]
        dj = 0
        for i = N0-N+1:N0
            @inbounds dj += isconnected(eg,eg.vertices[i],vj)
        end
        if dj < dmin
            dmin = dj
            vmin = vj
        end
    end
    return vmin
end

function mis2(graph::EliminateGraph, level::Int=0)
    if nv(graph) == 0
        return 0
    else
        v = minx(v->degree(graph,v), graph.vertices)
        return 1 + maximum(y->mis1(eliminate_neighborcover(graph,y), level+1), neighborcover(graph, v))
    end
end

@testset "mis1" begin
    @test minx(x->x^2, [2,-1, 3]) == -1

    graph = EliminateGraph([0 1 0 0 0;
                            1 0 1 1 0;
                            0 1 0 1 0;
                            0 1 1 0 1;
                            0 0 0 1 0])

    @test mis1(graph) == 3
    #@test mis2(graph) == 3
end

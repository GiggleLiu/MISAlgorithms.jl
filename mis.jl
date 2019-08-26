include("utils.jl")
include("EliminateGraph.jl")

function mis1(graph::EliminateGraph, level::Int=0)
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
end

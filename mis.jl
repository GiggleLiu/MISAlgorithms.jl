include("EliminateGraph.jl")

function mis1(graph::EliminateGraph, level::Int=0)
    if nv(graph) == 0
        return 0
    else
        v = minx(v->degree(graph,v), graph.vertices)
        return 1 + maximum(y->mis1(graph\neighborcover(graph, y), level+1), neighborcover(graph, v))
    end
end

function minx(f, vec)
    local xmin = vec[1]
    fmin = f(xmin)
    for j=2:length(vec)
        @inbounds x = vec[j]
        fmin_ = f(x)
        fmin_ < fmin && (fmin = fmin_; xmin=x)
    end
    return xmin
end

@test minx(x->x^2, [2,-1, 3]) == -1

graph = EliminateGraph([0 1 0 0 0;
                        1 0 1 1 0;
                        0 1 0 1 0;
                        0 1 1 0 1;
                        0 0 0 1 0])

@test mis1(graph) == 3

# mask-copy version, 1.5 second
using Random
Random.seed!(2)
eg = rand_egraph(50, 0.1)
using BenchmarkTools
@benchmark neighborcover($eg, 2) seconds=1


using Profile
Profile.clear()
Profile.init(delay=1e-3)
@profile neighborcover(eg, 2)

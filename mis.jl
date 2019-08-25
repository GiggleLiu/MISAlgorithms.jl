include("EliminateGraph.jl")

function mis1(graph::EliminateGraph, level::Int=0)
    if nv(graph) == 0
        return 0
    else
        vs = vertices(graph)
        v = minx(v->degree(graph,v,vs), vs)
        #v = minx(v->degree(graph,v), vs)
        #graph\neighborcover(graph, neighborcover(graph,v,vs)[1], vs)
        #return 1 + maximum(y->(@show level, y, vs; mis1(graph\neighborcover(graph, y, vs), level+1)), neighborcover(graph, v, vs))
        return 1 + maximum(y->mis1(graph\neighborcover(graph, y, vs), level+1), neighborcover(graph, v, vs))
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
eg = rand_egraph(60, 0.1)
@time mis1(eg)


using Profile
Profile.clear()
@profile for i=1:10000 mis1(eg) end

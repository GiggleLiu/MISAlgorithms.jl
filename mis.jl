include("EliminateGraph.jl")

function mis1(graph::EliminateGraph, level::Int=0)
    if nv(graph) == 0
        return 0
    else
        v = minx(v->iseliminated(graph, v) ? 999999 : degree(graph,v), 1:nv0(graph))
        res_ = maximum(NeighborCover(graph, v)) do y
            eliminate!(graph, NeighborCover(graph, y))
            ri = mis1(graph, level+1)
            recover!(graph, NeighborCover(graph, y))
            return ri
        end
        return 1 + res_
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

"""
    copyltu!(A::AbstractMatrix) -> AbstractMatrix

copy the lower triangular to upper triangular.
"""
function copyltu!(A::AbstractMatrix)
    m, n = size(A)
    for i=1:m
        A[i,i] = real(A[i,i])
        for j=i+1:n
            @inbounds A[i,j] = conj(A[j,i])
        end
    end
    A
end

function rand_egraph(nv::Int, density::Real)
    tbl = rand(nv, nv) .< density
    copyltu!(tbl)
    for i=1:nv
        tbl[i,i] = false
    end
    EliminateGraph(tbl)
end

using Random
Random.seed!(2)
eg = rand_egraph(60, 0.1)

# Bool Matrix: 1.9s
@time mis1(eg)

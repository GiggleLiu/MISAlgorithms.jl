include("EliminateGraph.jl")

function mis1(graph::EliminateGraph, level::Int=0)
    if nv(graph) == 0
        return 0
    else
        vs = vertices(graph)
        v = minx(v->degree(graph,v,vs), vs)
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
mis1(graph)

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

eg = rand_egraph(60, 0.1)

mis1(eg)

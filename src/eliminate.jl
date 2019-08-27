export eliminate, eliminate!, recover!

"""
    eliminate!(eg::EliminateGraph, vertices)

Eliminate vertices from a graph.
"""
function eliminate!(eg::EliminateGraph, vi::Int)
    N = nv0(eg)
    @inbounds iptr = eg.level == 0 ? 0 : eg.ptr[eg.level]
    for j in N-eg.nv+1:N
        @inbounds vj = eg.vertices[j]
        if vj==vi
            iptr += 1
            eg.level += 1
            @inbounds eg.ptr[eg.level] = iptr
            unsafe_swap!(eg.vertices, j, iptr)
            break
        end
    end
    eg.nv -= 1
    return eg
end

function eliminate!(eg::EliminateGraph, vs)
    N = nv0(eg)
    @inbounds iptr = eg.level == 0 ? 0 : eg.ptr[eg.level]
    for vi in vs
        for j in N-eg.nv+1:N
            @inbounds vj = eg.vertices[j]
            if vj==vi
                iptr += 1
                unsafe_swap!(eg.vertices, j, iptr)
                break
            end
        end
    end
    eg.level += 1
    eg.nv -= length(vs)
    @inbounds eg.ptr[eg.level] = iptr
    return eg
end

function eliminate!(eg::EliminateGraph, nc::NeighborCover)
    N = nv0(eg)
    @inbounds iptr = eg.level == 0 ? 0 : eg.ptr[eg.level]
    iptr0 = iptr
    for i in N-eg.nv+1:N
        @inbounds vi = eg.vertices[i]
        (vi == nc.i || isconnected(eg, vi, nc.i)) || continue
        for j in N-eg.nv+1:N
            @inbounds vj = eg.vertices[j]
            if vj==vi
                iptr += 1
                eg.nv -= 1
                unsafe_swap!(eg.vertices, j, iptr)
                break
            end
        end
    end
    eg.level += 1
    @inbounds eg.ptr[eg.level] = iptr
    return eg
end

"""
    eliminate([func], eg::EliminateGraph, vertices)
    eg \\ vertices

Eliminate vertices from a graph, return the value of `func(eliminated_graph)` if `func` provided.
"""
eliminate(eg, vertices) = eliminate!(copy(eg), vertices)

@inline function eliminate(func, eg::EliminateGraph, vi)
    eliminate!(eg, vi)
    res = func(eg)
    recover!(eg)
    return res
end

Base.:\(eg::EliminateGraph, vertices) = eliminate(eg, vertices)

"""restore eliminated vertices for a level (one call of elimintion function)."""
function recover!(eg::EliminateGraph)
    @inbounds eg.nv += eg.ptr[eg.level] - (eg.level==1 ? 0 : eg.ptr[eg.level-1])
    eg.level -= 1
    eg
end

export neighborcover, neighborcover_mapreduce, mirrorcover, mirrors, isclique

"""
    neighborcover(eg::EliminateGraph, vi::Int) -> EliminateGraph

Return `vi` and its neighbors.
"""
@inline function neighborcover(eg::EliminateGraph, vi::Int)
    return filter(vj->isconnected(eg,vj,vi) || vi==vj, vertices(eg))
end

function neighborcover_mapreduce(func, red, eg::EliminateGraph, vi::Int)
    pre = func(vi)
    for vj in vertices(eg)
        if isconnected(eg,vj,vi)
            # find a neighbor cover
            pre = red(pre, func(vj))
        end
    end
    return pre
end

function isclique(eg::EliminateGraph, vs)
    for vj in vs
        for vi in vs
            vi!=vj && !isconnected(eg, vi, vj) && return false
        end
    end
    return true
end

isclique(eg::EliminateGraph) = isclique(eg, vertices(eg))

function mirrors(eg::EliminateGraph, vi::Int, nbs2)
    nc = neighborcover(eg, vi)
    filter(v->isclique(eg, setdiff(nc,v)), nbs2)
end

function mirrors(eg::EliminateGraph, vi::Int)
    mirrors(eg, vi, neighbors2(eg, vi))
end

function mirrorcover(eg::EliminateGraph, vi::Int, nbs2)
    nc = neighborcover(eg, vi)
    push!(filter(v->isclique(eg, setdiff(nc,neighbors(eg,v))), nbs2), vi)
end

function mirrorcover(eg::EliminateGraph, vi::Int)
    mirrorcover(eg, vi, neighbors2(eg, vi))
end

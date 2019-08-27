export neighbors, degrees, degree, mindegree_vertex, maxdegree_vertex, minmaxdegree_vertex

function neighbors(eg::EliminateGraph, i::Int)
    return filter(j->isconnected(eg,j,i), vertices(eg))
end

degrees(eg::EliminateGraph) = degree.(Ref(eg), vertices(eg))

@inline function degree(eg::EliminateGraph, vi::Int)
    sum(vj->isconnected(eg,vi,vj), vertices(eg))
end

function mindegree_vertex(eg::EliminateGraph)
    N = nv(eg)
    dmin = 999999
    vmin = 0
    for vj in vertices(eg)
        dj = 0
        for vi in vertices(eg)
            dj += isconnected(eg,vi,vj)
        end
        if dj < dmin
            dmin = dj
            vmin = vj
        end
    end
    return vmin, dmin
end

function minmaxdegree_vertex(eg::EliminateGraph)
    N = nv(eg)
    dmin = 999999
    dmax = 0
    vmin = 0
    vmax = 0
    for vj in vertices(eg)
        dj = 0
        for vi in vertices(eg)
            dj += isconnected(eg,vi,vj)
        end
        if dj < dmin
            dmin = dj
            vmin = vj
        end
        if dj > dmax
            dmax = dj
            vmax = vj
        end
    end
    return vmin, vmax, dmin, dmax
end

function maxdegree_vertex(eg::EliminateGraph)
    N = nv(eg)
    dmax = 0
    vmax = 0
    for vj in vertices(eg)
        dj = 0
        for vi in vertices(eg)
            dj += isconnected(eg,vi,vj)
        end
        if dj > dmax
            dmax = dj
            vmax = vj
        end
    end
    return vmax, dmax
end

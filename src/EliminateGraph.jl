export EliminateGraph, rand_eg, refresh
export nv0, nv, isconnected

"""
    EliminateGraph
    EliminateGraph(tbl::AbstractMatrix) -> EliminateGraph

A graph type for algorithms that involve node elimination.
With this type, vertex elimination and recover do not allocate.
`tbl` in the constructor is a boolean table for connection.
"""
mutable struct EliminateGraph
    tbl::Matrix{Bool}
    vertices::Vector{Int}
    ptr::Vector{Int}
    level::Int
    nv::Int
end

EliminateGraph(tbl::AbstractMatrix) = EliminateGraph(Matrix{Bool}(tbl))
function EliminateGraph(tbl::Matrix{Bool})
    N = size(tbl, 1)
    vertices = collect(1:N)
    ptr = zeros(Int,N)
    EliminateGraph(tbl, vertices, ptr, 0, N)
end

Base.copy(eg::EliminateGraph) = EliminateGraph(eg.tbl, eg.vertices |> copy, eg.ptr|>copy, eg.level, eg.nv)

"""initial size of a `EliminateGraph`."""
nv0(eg::EliminateGraph) = size(eg.tbl, 1)
"""current size of a `EliminateGraph`."""
nv(eg::EliminateGraph) = eg.nv

"""undo elimination for a `EliminateGraph`."""
refresh(eg::EliminateGraph) = EliminateGraph(eg.tbl)

"""
    isconnected(eg::EliminateGraph, i::Int, j::Int) -> Bool

Return true if `i`, `j` are connected in `eg`.
Note: This function does not check `i`, `j` out of bound error!
"""
isconnected(eg::EliminateGraph, i::Int, j::Int) = @inbounds eg.tbl[i,j]

function Base.show(io::IO, eg::EliminateGraph)
    N = nv0(eg)
    println(io, "EliminateGraph")
    vs = vertices(eg)
    for i=1:N
        for j=1:N
            print(io, "  ", (i in vs && j in vs) ? Int(eg.tbl[i,j]) : "⋅")
        end
        println(io)
    end
end

"""
    rand_eg(nv::Int, density::Real) -> EliminateGraph

Generate a random `EliminateGraph`.
"""
function rand_eg(nv::Int, density::Real)
    tbl = rand(nv, nv) .< density
    copyltu!(tbl)
    for i=1:nv
        tbl[i,i] = false
    end
    EliminateGraph(tbl)
end

export Vertices, vertices, EliminatedVertices, eliminated_vertices, iseliminated

"""
    Vertices{GT}<:AbstractArray{Int,1}

Vertex enumerator for a graph.
"""
struct Vertices{GT}<:AbstractArray{Int,1}
    eg::GT
end
"""vertices of a graph."""
vertices(eg::EliminateGraph) = Vertices(eg)
Base.length(vs::Vertices) = vs.eg.nv
Base.getindex(vs::Vertices, i::Int) = vs.eg.vertices[end-vs.eg.nv+i]

"""
    EliminatedVertices{GT}<:AbstractArray{Int,1}

Eliminated vertex enumerator for a graph.
"""
struct EliminatedVertices{GT}<:AbstractArray{Int,1}
    eg::GT
end
Base.length(vs::EliminatedVertices) = nv0(vs.eg)-vs.eg.nv
Base.getindex(vs::EliminatedVertices, i::Int) = vs.eg.vertices[i]

"""eliminated vertices of a `EliminateGraph`."""
eliminated_vertices(eg::EliminateGraph) = EliminatedVertices(eg)

"""
    iseliminated(eg::EliminateGraph, i::Int) -> Bool

Return true if a vertex of a `EliminateGraph` is eliminated.
"""
iseliminated(eg::EliminateGraph, i::Int) = i in eliminated_vertices(eg)

for V in [:Vertices, :EliminatedVertices]
    @eval Base.size(vs::$V) = (length(vs),)
    @eval Base.size(vs::$V, i::Int) = i==1 ? length(vs) : 1
    @eval Base.IteratorEltype(::Type{$V}) = Base.HasEltype()
    @eval Base.IteratorSize(::Type{$V}) = Base.HasLength()
    @eval Base.eltype(vs::$V) = Int
    @eval Base.eltype(vs::Type{$V}) = Int

    VI = V==:Vertices ? :(eg.vertices[end-eg.nv+state]) : :(eg.vertices[state])
    @eval function Base.iterate(vs::$V, state=1)
        eg = vs.eg
        if state > length(vs)
            return nothing
        else
            @inbounds vi = $VI
            return vi, state+1
        end
    end
end

export neighborcover, neighborcover_mapreduce, mirrorcover

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

function mirrorcover(eg::EliminateGraph, vi::Int)
end

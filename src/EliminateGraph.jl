export EliminateGraph, rand_egraph, refresh
export nv0, nv, isconnected

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

nv0(eg::EliminateGraph) = size(eg.tbl, 1)
Base.copy(eg::EliminateGraph) = EliminateGraph(eg.tbl, eg.vertices |> copy, eg.ptr|>copy, eg.level, eg.nv)
refresh(eg::EliminateGraph) = EliminateGraph(eg.tbl)

isconnected(eg::EliminateGraph, i::Int, j::Int) = @inbounds eg.tbl[i,j]
nv(eg::EliminateGraph) = eg.nv

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

function rand_egraph(nv::Int, density::Real)
    tbl = rand(nv, nv) .< density
    copyltu!(tbl)
    for i=1:nv
        tbl[i,i] = false
    end
    EliminateGraph(tbl)
end

export Vertices, vertices, EliminatedVertices, eliminated_vertices, iseliminated

struct Vertices{GT}<:AbstractArray{Int,1}
    eg::GT
end
vertices(eg::EliminateGraph) = Vertices(eg)
Base.length(vs::Vertices) = vs.eg.nv
Base.getindex(vs::Vertices, i::Int) = vs.eg.vertices[end-vs.eg.nv+i]

struct EliminatedVertices{GT}<:AbstractArray{Int,1}
    eg::GT
end
Base.length(vs::EliminatedVertices) = nv0(vs.eg)-vs.eg.nv
Base.getindex(vs::EliminatedVertices, i::Int) = vs.eg.vertices[i]
eliminated_vertices(eg::EliminateGraph) = EliminatedVertices(eg)
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

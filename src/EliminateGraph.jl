export EliminateGraph, subgraph, refresh, check_validity
export nv0, nv, ne0, ne
export isconnected, unsafe_connect!, find_cluster

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
function EliminateGraph(nv::Int, pairs::Vector{Pair{Int,Int}})
    tbl = zeros(Bool, nv, nv)
    for (vi, vj) in pairs
        tbl[vi,vj] = true
    end
    return EliminateGraph(tbl .| tbl')
end

Base.copy(eg::EliminateGraph) = EliminateGraph(eg.tbl, eg.vertices |> copy, eg.ptr|>copy, eg.level, eg.nv)
Base.:(==)(eg1::EliminateGraph, eg2::EliminateGraph) = eg1.tbl == eg2.tbl && Set(vertices(eg1)) == Set(vertices(eg2))

"""initial size of a `EliminateGraph`."""
nv0(eg::EliminateGraph) = size(eg.tbl, 1)
"""initial number of edges in a `EliminateGraph`."""
ne0(eg::EliminateGraph) = sum(eg.tbl) ÷ 2
"""current size of a `EliminateGraph`."""
nv(eg::EliminateGraph) = eg.nv
"""current number of edges in a `EliminateGraph`."""
function ne(eg::EliminateGraph)
    res = 0
    for vj in vertices(eg)
        for vi in vertices(eg)
            vi > vj && (res += isconnected(eg, vi, vj))
        end
    end
    return res
end

"""
A eliminate graph is valid.
"""
function check_validity(eg::EliminateGraph)
    mes = ""
    eg.tbl' == eg.tbl || (mes *= "connection table not symmetric\n")
    !any(i->eg.tbl[i,i], 1:nv0(eg)) || (mes *= "diagonal part of table non-zero\n")
    Set(eg.vertices) == Set(1:nv0(eg)) || (mes *= "vertices list incorrect\n")
    0 <= eg.nv <= nv0(eg) || (mes *= "number of vertices not in range\n")
    0 <= eg.level <= nv0(eg) || (mes *= "level not in range\n")
    if eg.level > 0
        eg.ptr[eg.level] == nv0(eg)-eg.nv || (mes *= "ptr points to incorrect position\n")
        nv0(eg) >= eg.ptr[1] > 0 || (mes *= "ptr not in range\n")
        all(diff(eg.ptr[1:eg.level]) .> 0) || (mes *= "non-increasing ptr")
    end
    if mes == ""
        return true
    else
        println(mes)
        return false
    end
end

"""undo elimination for a `EliminateGraph`."""
refresh(eg::EliminateGraph) = EliminateGraph(eg.tbl)

"""
    isconnected(eg::EliminateGraph, vi::Int, vj::Int) -> Bool

Return true if `vi`, `vj` are connected in `eg`.
Note: This function does not check `vi`, `vj` out of bound error!
"""
isconnected(eg::EliminateGraph, vi::Int, vj::Int) = @inbounds eg.tbl[vi,vj]

getid(eg::EliminateGraph, vi::Int) = findfirst(==(vi), vertices(eg))

"""
    find_cluster(eg::EliminateGraph, vi::Int) -> Vector{Int}

Find the cluster connected to `vi` in `eg`.
"""
function find_cluster(eg::EliminateGraph, vi::Int)
    vset = zeros(Bool, nv0(eg))
    find_cluster!(eg, vi, vset)
    return findall(!iszero, vset)
end

function find_cluster!(eg::EliminateGraph, vi::Int, vset::Vector{Bool})
    vset[vi] && return

    vset[vi] = true
    for vj in neighbors(eg, vi)
        find_cluster!(eg, vj, vset)
    end
end

function isconnected(eg::EliminateGraph)
    eg.nv <=1 && return true
    length(find_cluster(eg, vertices(eg)[1])) == nv(eg)
end

subgraph(eg::EliminateGraph, vs) = EliminateGraph(eg.tbl[vs, vs])

"""
    unsafe_connect!(eg::EliminateGraph, vi::Int, vj::Int) -> EliminateGraph

connect two vertices.
"""
function unsafe_connect!(eg::EliminateGraph, vi::Int, vj::Int)
    @inbounds eg.tbl[vi, vj] = true
    @inbounds eg.tbl[vj, vi] = true
    eg
end

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

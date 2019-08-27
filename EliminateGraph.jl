abstract type AbstractVertexSet end
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

vertices(eg::EliminateGraph) = Vertices(eg)
eliminated_vertices(eg::EliminateGraph) = EliminatedVertices(eg)
iseliminated(eg::EliminateGraph, i::Int) = i in eliminated_vertices(eg)
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

struct NeighborCover <: VertexSet
    i::Int
end

@inline function unsafe_swap!(v::Vector, i::Int, j::Int)
    @inbounds temp = v[i]
    @inbounds v[i] = v[j]
    @inbounds v[j] = temp
end

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

eliminate(eg, vertices) = eliminate!(copy(eg), vertices)

@inline function eliminate(func, eg::EliminateGraph, vi)
    eliminate!(eg, vi)
    res = func(eg)
    recover!(eg)
    return res
end

Base.:\(eg::EliminateGraph, vertices) = eliminate(eg, vertices)
degrees(eg::EliminateGraph) = degrees(eg, vertices(eg))
function degrees(eg::EliminateGraph)
    map(vi->degree(eg, vi), vertices(eg))
end

@inline function degree(eg::EliminateGraph, vi::Int)
    sum(vj->isconnected(eg,vi,vj), vertices(eg))
end

function recover!(eg::EliminateGraph)
    @inbounds eg.nv += eg.ptr[eg.level] - (eg.level==1 ? 0 : eg.ptr[eg.level-1])
    eg.level -= 1
    eg
end

function neighbors(eg::EliminateGraph, i::Int)
    return filter(j->!iseliminated(eg,j) && isconnected(eg,j,i), 1:nv0(eg))
end
function neighbors(eg::EliminateGraph, i::Int)
    return filter(j->isconnected(eg,j,i), vertices(eg))
end

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

struct Vertices{GT}<:AbstractArray{Int,1}
    eg::GT
end

Base.length(vs::Vertices) = vs.eg.nv
Base.size(vs::Vertices) = (length(vs),)
Base.size(vs::Vertices, i::Int) = i==1 ? length(vs) : throw(DimensionMismatch(""))
Base.eltype(vs::Vertices) = Int
Base.IteratorEltype(::Type{Vertices}) = Base.HasEltype()
Base.IteratorSize(::Type{Vertices}) = Base.HasLength()
Base.eltype(vs::Type{Vertices}) = Int

function Base.iterate(vs::Vertices, state=1)
    eg = vs.eg
    if state > eg.nv
        return nothing
    else
        @inbounds vi = eg.vertices[end-eg.nv+state]
        return vi, state+1
    end
end

Base.getindex(vs::Vertices, i::Int) = vs.eg.vertices[end-vs.eg.nv+i]

struct EliminatedVertices{GT}<:AbstractArray{Int,1}
    eg::GT
end

Base.length(vs::EliminatedVertices) = nv0(vs.eg)-vs.eg.nv
Base.size(vs::EliminatedVertices) = (length(vs),)
Base.size(vs::EliminatedVertices, i::Int) = i==1 ? length(vs) : throw(DimensionMismatch(""))
Base.eltype(vs::EliminatedVertices) = Int

function Base.iterate(vs::EliminatedVertices, state=1)
    eg = vs.eg
    if state > length(vs)
        return nothing
    else
        @inbounds vi = eg.vertices[state]
        return vi, state+1
    end
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

struct Vertex<:AbstractVertexSet
    i::Int
end

struct UnionOf{TA, TB}<:AbstractVertexSet
    A::TA
    B::TB
end
∪(A::Union{Vertex, NeighborCover, MirrorCover}, B::)

function mirrors(eg, v)
    Vertex()
end

using Test
@testset "eliminate and recover" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @show eg
    @test neighbors(eg, 1) == [2]
    @test neighborcover(eg, 1) == [1,2]
    @test nv(eg) == 4
    eliminate!(eg, 2)
    eliminate!(eg, 3)
    @test nv(eg) == 2
    recover!(eg)
    @test nv(eg) == 3
    eg2 = eg \ 3
    @test nv(eg) == 3
    @test nv(eg2) == 2
    @test vertices(eg2) == [1,4]
    @test neighbors(eg2, 1) == []
    @test degrees(eg) == [1, 0, 1]
    @test degrees(eg2) == [0, 0]
    eg3 = eg2 \ (1,4)
    @test nv(eg3) == 0

    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @show eg
    res = eliminate(eg, 3) do eg
        nv(eg)
    end # do
    @test res == 3

    @show eg
    eg4 = eg \ NeighborCover(1)
    @show eg4
    @test vertices(eg4) == [3,4]

    res = eliminate(eg, (3,4)) do eg
        eliminate(eg, (1,2)) do eg
            nv(eg)
        end
    end
    @test res == 0
    @test nv(eg) == 4
    @test eg.level == 0

    vs = Vertices(eg)
    @test [vs...] == vertices(eg)
    @test vs[3] == vertices(eg)[3]
    vmin, dmin = mindegree_vertex(eg)
    @test degree(eg, vmin) == minimum(degrees(eg)) == dmin
    vmax, dmax = maxdegree_vertex(eg)
    @test degree(eg, vmax) == maximum(degrees(eg)) == dmax
    vmin, vmax, dmin, dmax = minmaxdegree_vertex(eg)
    @test degree(eg,vmin) == minimum(degrees(eg)) == dmin
    @test degree(eg,vmax) == maximum(degrees(eg)) == dmax
end

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

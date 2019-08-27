export AbstractVertexSet, NeighborCover, MirrorCover, Vertex, UnionOf, Vertex

# Vertex Set Notations
abstract type AbstractVertexSet end

struct NeighborCover <: AbstractVertexSet
    i::Int
end

struct MirrorCover <: AbstractVertexSet
    i::Int
end

struct Vertex<:AbstractVertexSet
    i::Int
end

struct UnionOf{TA, TB}<:AbstractVertexSet
    A::TA
    B::TB
end
Base.:∪(A::AbstractVertexSet, B::AbstractVertexSet) = UnionOf(A, B)

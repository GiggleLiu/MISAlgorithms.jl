export AbstractVertexSet, NeighborCover, MirrorCover, Vertex, UnionOf, Vertex

# Vertex Set Notations
abstract type AbstractVertexSet end

"""
    NeighborCover <: AbstractVertexSet

Neighbors and itself.
"""
struct NeighborCover <: AbstractVertexSet
    i::Int
end

"""
    MirrorCover <: AbstractVertexSet

Mirrors and itself. A vertex `w` is a mirror of `v` if `w ∈ N²(v)` and `N[v]\\N(w)` is a clique.
"""
struct MirrorCover <: AbstractVertexSet
    i::Int
end

struct Vertex<:AbstractVertexSet
    i::Int
end

"""
    UnionOf{TA, TB}<:AbstractVertexSet

A union of ...
"""
struct UnionOf{TA, TB}<:AbstractVertexSet
    A::TA
    B::TB
end
Base.:∪(A::AbstractVertexSet, B::AbstractVertexSet) = UnionOf(A, B)

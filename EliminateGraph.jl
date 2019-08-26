struct EliminateGraph
    tbl::Matrix{Bool}
    vertices::Vector{Int}
end

nv0(eg::EliminateGraph) = size(eg.tbl, 1)

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

EliminateGraph(tbl::AbstractMatrix) = EliminateGraph(Matrix{Bool}(tbl))
EliminateGraph(tbl::Matrix{Bool}) = EliminateGraph(tbl, collect(1:size(tbl,1)))

function eliminate!(eg::EliminateGraph, vi::Int)
    for i in 1:nv(eg)
        if vi == eg.vertices[i]
            deleteat!(eg.vertices,i)
            return eg
        end
    end
end

function eliminate!(eg::EliminateGraph, vertices)
    for vi in vertices
        eliminate!(eg, vi)
    end
    eg
end

Base.copy(eg::EliminateGraph) = EliminateGraph(eg.tbl, eg.vertices |> copy)

@inline function eliminate(eg::EliminateGraph, vertices)
    vs = Int[]
    @inbounds rmv = vertices[1]
    rmk = 1
    for iv in eg.vertices
        if iv == rmv
            rmk+=1
            rmk <= length(vertices) && (rmv = @inbounds vertices[rmk])
        else
            push!(vs, iv)
        end
    end
    EliminateGraph(eg.tbl, vs)
end

@inline function eliminate(eg::EliminateGraph, vertices)
    vs = Int[]
    @inbounds rmv = vertices[1]
    rmk = 1
    for iv in eg.vertices
        if iv == rmv
            rmk+=1
            rmk <= length(vertices) && (rmv = @inbounds vertices[rmk])
        else
            push!(vs, iv)
        end
    end
    EliminateGraph(eg.tbl, vs)
end

Base.:\(eg::EliminateGraph, vertices) = eliminate(eg, vertices)
degrees(eg::EliminateGraph) = degrees(eg, vertices(eg))
function degrees(eg::EliminateGraph)
    map(vi->degree(eg, vi), vertices(eg))
end

@inline function degree(eg::EliminateGraph, vi::Int)
    sum(vj->isconnected(eg,vi,vj), vertices(eg))
end

function recover!(eg::EliminateGraph, vi::Int)
    push!(eg.vertices, vi)
    eg
end

vertices(eg::EliminateGraph) = eg.vertices
iseliminated(eg::EliminateGraph, i::Int) = !(i in vertices(eg))
isconnected(eg::EliminateGraph, i::Int, j::Int) = @inbounds eg.tbl[i,j]
nv(eg::EliminateGraph) = length(vertices(eg))
function neighbors(eg::EliminateGraph, i::Int)
    return filter(j->!iseliminated(eg,j) && isconnected(eg,j,i), 1:nv0(eg))
end
function neighbors(eg::EliminateGraph, i::Int)
    return filter(j->isconnected(eg,j,i), vertices(eg))
end

@inline function neighborcover(eg::EliminateGraph, vi::Int)
    return filter(vj->isconnected(eg,vj,vi) || vi==vj, vertices(eg))
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
    recover!(eg, 3)
    @test nv(eg) == 3
    eg2 = eg \ 3
    @test nv(eg) == 3
    @test nv(eg2) == 2
    @test vertices(eg2) == [1,4]
    @test neighbors(eg2, 1) == []
    @test degrees(eg) == [0, 1, 1]
    @test degrees(eg2) == [0, 0]
    eg3 = eg2 \ (1,4)
    @test nv(eg3) == 0
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

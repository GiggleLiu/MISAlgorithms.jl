struct EliminateGraph
    tbl::BitMatrix
    mask::BitVector
end

nv0(eg::EliminateGraph) = size(eg.tbl, 1)

function Base.show(io::IO, eg::EliminateGraph)
    N = nv0(eg)
    println(io, "EliminateGraph")
    for i=1:N
        for j=1:N
            print(io, "  ", (eg.mask[i] && eg.mask[j]) ? Int(eg.tbl[i,j]) : "⋅")
        end
        println(io)
    end
end

EliminateGraph(tbl::Matrix) = EliminateGraph(BitMatrix(tbl))
EliminateGraph(tbl::BitMatrix) = EliminateGraph(tbl, BitVector(ones(Bool,size(tbl,1))))

function eliminate!(eg::EliminateGraph, i::Int)
    eg.mask[i] = false
    eg
end

function eliminate!(eg::EliminateGraph, vertices)
    for i in vertices
        eg.mask[i] = false
    end
    eg
end

Base.copy(eg::EliminateGraph) = EliminateGraph(eg.tbl, eg.mask |> copy)
eliminate(eg::EliminateGraph, vertices) = eliminate!(copy(eg), vertices)
Base.:\(eg::EliminateGraph, vertices) = eliminate(eg, vertices)
degrees(eg::EliminateGraph) = degrees(eg, vertices(eg))
function degrees(eg::EliminateGraph, vertices)
    map(vi->degree(eg, vi, vertices), vertices)
end

@inline function degree(eg::EliminateGraph, vi::Int, vertices)
    sum(vj->isconnected(eg,vi,vj), vertices)
end

function recover!(eg::EliminateGraph, i::Int)
    eg.mask[i] = true
    eg
end

vertices(eg::EliminateGraph) = findall(eg.mask)
iseliminated(eg::EliminateGraph, i::Int) = !eg.mask[i]
isconnected(eg::EliminateGraph, i::Int, j::Int) = eg.tbl[i,j]
nv(eg::EliminateGraph) = sum(eg.mask)
function neighbors(eg::EliminateGraph, i::Int)
    return filter(j->!iseliminated(eg,j) && isconnected(eg,j,i), 1:nv0(eg))
end
function neighbors(eg::EliminateGraph, i::Int, vertices)
    return filter(j->isconnected(eg,j,i), vertices)
end

function neighborcover(eg::EliminateGraph, i::Int, vertices)
    return filter(j->(i==j || isconnected(eg,j,i)), vertices)
end

using Test
@testset "eliminate and recover" begin
    tbl = Bool[false true true false; true false true true; false true false true; true true true false]
    eg = EliminateGraph(tbl .& tbl')
    @show eg
    @test neighbors(eg, 1) == [2]
    @test neighborcover(eg, 1, vertices(eg)) == [1,2]
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

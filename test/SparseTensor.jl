using BitBasis, SparseArrays

# BitStr patch
function Base.sub_with_overflow(x::T, y::T) where T<:BitStr
    res, sig = Base.sub_with_overflow(buffer(x), buffer(y))
    return T(res), sig
end

function Base.add_with_overflow(x::T, y::T) where T<:BitStr
    res, sig = Base.add_with_overflow(buffer(x), buffer(y))
    return T(res), sig
end

function bpermute(b::T, order) where T<:Integer
    res = zero(b)
    for i = 1:length(order)
        @inbounds bi = order[i]
        res |= (b & bmask(T,bi)) >> (bi-i)
    end
    return res
end

struct BinarySparseTensor{Tv,Ti<:Integer,N} <: AbstractSparseArray{Tv, Ti, N}
   data:: SparseVector{Tv, Ti}
end

as_index(x::Integer) = x
as_index(x::BitStr) = buffer(x)+1
@inline function as_index(x::NTuple{N,<:Integer}) where N
    res = 1
    for i=1:N
        @inbounds res += (x[i]-1)<<(i-1)
    end
    return res
end

function bst(data::SparseVector{T,Ti}) where {T,Ti}
    N = log2dim1(data)
    length(data) != 1<<N && throw(ArgumentError("data length should be 2^x, got $(length(data))"))
    BinarySparseTensor{T,Ti,N}(data)
end

function BinarySparseTensor(A::AbstractArray)
    bst(SparseVector(vec(A)))
end
Base.getindex(t::BinarySparseTensor{T,Ti,N}, index::BitStr{N,Ti}) where {T,Ti,N} = @inbounds t.data[as_index(index)]
Base.getindex(t::BinarySparseTensor, index::Int...) = t.data[as_index(index)]
function Base.size(t::BinarySparseTensor{T,Ti,N}, i::Int) where {T,Ti,N}
    i<=0 && throw(BoundsError(size(t), i))
    i<=N ? 2 : 1
end
Base.size(t::BinarySparseTensor{T,Ti,N}) where {T,Ti,N} = (fill(2, N)...,)

function Base.setindex!(t::BinarySparseTensor{T,Ti,N}, val, index::BitStr{N,Ti}) where {T,Ti,N}
    @inbounds t.data[as_index(index)] = val
end

SparseArrays.nnz(t::BinarySparseTensor) = nnz(t.data)
SparseArrays.findnz(t::BinarySparseTensor{Tv,Ti,N}) where {Tv,Ti,N} = findnz(BitStr{N,Ti}, t)
function SparseArrays.findnz(::Type{T}, t::BinarySparseTensor) where T
    convert.(T, t.data.nzind.-1), t.data.nzval
end

Base.show(io::IO, ::MIME"text/plain", t::BinarySparseTensor) = Base.show(io, t)
function Base.show(io::IO, t::BinarySparseTensor{T,Ti,N}) where {T,Ti,N}
    NNZ = length(t.data.nzind)
    println(io, "$(summary(t)) with $(nnz(t)) stored entries:")
    for (i, (nzi, nzv)) in enumerate(zip(findnz(t)...))
        print(io, "  $nzi = $nzv")
        i != NNZ && println(io)
    end
end

bst_zeros(::Type{Tv}, ::Type{Ti}, N::Int) where {Tv,Ti} = BinarySparseTensor{Tv,Ti,N}(SparseVector(1<<N, Ti[], Tv[]))
bst_zeros(::Type{Tv}, N::Int) where {Tv} = bst_zeros(Tv, Int64, N)

Base.zero(t::BinarySparseTensor) = bst(SparseVector(t.data.n, t.data.nzind, zero(t.data.nzval)))

function Base.permutedims!(dest::BinarySparseTensor{Tv,Ti,N}, src::BinarySparseTensor{Tv,Ti,N}, dims::NTuple{N,Int}) where {Tv,Ti,N}
    nzind, nzval = findnz(src)
    nzind .= bpermute.(nzind, Ref(dims))
    order = sortperm(nzind)
    @inbounds dest.data.nzind .= buffer.(nzind[order]) .+ 1
    @inbounds dest.data.nzval .= nzval[order]
    return dest
end

function Base.permutedims(src::BinarySparseTensor{Tv,Ti,N}, dims::NTuple{N,Int}) where {Tv,Ti,N}
    dest = zero(src)
    Base.permutedims!(dest, src, dims)
end

#function Base.reshape(a::BinarySparseTensor, dims::NTuple{N,Int}) where N
    #throw(DimensionMismatch("reshaping a BinarySparseTensor is not allowed."))
#end

using Test
@testset "constructor" begin
    sv = SparseVector([1,0,0,1,1,0,0,0])
    t = bst(sv)
    @test ndims(t) == 3
    @test size(t) == (2,2,2)
    @test nnz(t) == 3
    t2 = BinarySparseTensor(Array(t))
    @test t2 == t
    @test nnz(t2) == 3
    @test_throws BoundsError size(t,-1)
    @test size(t,3) == 2
    @test size(t,4) == 1
    @test size(t) == (2,2,2)
    cis = CartesianIndices((2,2,2))
    @test t[bit"000"] == 1
    @test t[bit"001"] == t[2,1,1] == 0 == t[cis[2,1,1]]
    @test t[bit"100"] == t[1,1,2] == 1 == t[cis[1,1,2]]
    @test permutedims(t, (2,1,3)) isa BinarySparseTensor
    @test permutedims(t, (2,1,3)) == permutedims(Array(t), (2,1,3))

    m = zeros(Int, 2,2,2); m[[1,4,5]] .= 1
    @test Array(t) == m
    t[bit"001"] = 8
    @test t[bit"001"] == 8
    sv = SparseVector([1,0,0,1,1,0,0,0,1])
    @test findnz(t) == ([bit"000", bit"001", bit"011", bit"100"], [1,8,1,1])
    @test_throws ArgumentError bst(sv)

    println(t)
    @test zero(t) == zeros(size(t))

    t = bst_zeros(Float64, 5)
    @test size(t) == (2,2,2,2,2)
    @test eltype(t) == Float64
    #@test_throws DimensionMismatch reshape(t, (2,2,2,2,2))
end

@testset "sort and permute" begin
    @test bpermute(0b1100, [2,4,3,1]) === 0b0110
    @test bpermute(0b1100, [4,2,3,1]) === 0b0101
    @test bpermute(bit"1100", [2,4,3,1]) === bit"0110"
    @test bpermute(bit"1100", [4,2,3,1]) === bit"0101"
    @test sort(BitStr64{5}[1,8,4,2,9]) == sort([1,8,4,2,9])
end

function bstrand(ndim::Int, density::Real)
    bst(sprand(1<<ndim, density))
end

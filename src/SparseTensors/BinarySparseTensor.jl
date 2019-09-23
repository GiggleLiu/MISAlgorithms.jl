using BitBasis, SparseArrays
using Base.Cartesian
using OMEinsum

export bst_zeros, bstrand, BinarySparseTensor, sparse_contract, bst, bpermute

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
SparseArrays.nonzeroinds(t::BinarySparseTensor{Tv,Ti,N}) where {Tv, Ti, N} = convert.(BitStr{N,Ti}, t.data.nzind.-1)
SparseArrays.nonzeros(t::BinarySparseTensor{Tv,Ti,N}) where {Tv, Ti, N} = t.data.nzval
Base.Array(t::BinarySparseTensor{Tv,Ti,1}) where {Tv,Ti} = Base.Array(t.data)

Base.show(io::IO, ::MIME"text/plain", t::BinarySparseTensor) = Base.show(io, t)
function Base.show(io::IOContext, t::BinarySparseTensor{T,Ti,1}) where {T,Ti,N}
    invoke(show, Tuple{IO,BinarySparseTensor}, io, t)
end
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
Base.copy(t::BinarySparseTensor) = bst(SparseVector(t.data.n, copy(t.data.nzind), copy(t.data.nzval)))

function Base.permutedims!(dest::BinarySparseTensor{Tv,Ti,N}, src::BinarySparseTensor{Tv,Ti,N}, dims::NTuple{N,Int}) where {Tv,Ti,N}
    nzind, nzval = findnz(src)
    nzind .= bpermute.(nzind, Ref(dims))
    order = sortperm(nzind)
    @inbounds dest.data.nzind .= buffer.(nzind[order]) .+ 1
    @inbounds dest.data.nzval .= nzval[order]
    return dest
end

Base.permutedims!(dest::BinarySparseTensor{Tv,Ti,N}, src::BinarySparseTensor{Tv,Ti,N}, dims) where {Tv,Ti,N} = permutedims!(dest, src, (dims...,))
function Base.permutedims(src::BinarySparseTensor{Tv,Ti,N}, dims) where {Tv,Ti,N}
    dest = copy(src)
    Base.permutedims!(dest, src, dims)
end

function bstrand(ndim::Int, density::Real)
    bst(sprand(1<<ndim, density))
end

include("patches.jl")
include("batched_gemm.jl")
include("index_reduction.jl")

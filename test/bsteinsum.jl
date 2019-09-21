using OMEinsum
include("SparseTensor.jl")

OMEinsum.asarray(x::Number, arr::BinarySparseTensor) where T = BinarySparseTensor(asarray(x))

function get_output_array(xs::NTuple{N, BinarySparseTensor{Tv,Ti,M} where {Tv,M}}, size) where {N,Ti}
    out = bst_zeros(promote_type(map(eltype,xs)...), Ti, length(size))
end

function loop_einsum!(code::EinCode{ixs, iy},
                xs::NTuple{N, BinarySparseTensor},
                y::BinarySparseTensor, size_dict) where {N, ixs, iy}
    NO = length(tunique(iy))
    A = einarray(code, xs, size_dict)
    y = reshape(y, fill(1, ndims(A)-NO)...,size(y)...)
    dropdims(Base._mapreducedim!(x->x, +, y, A), dims=(1:ndims(A)-NO...,))
end

@inline Base.getindex(A::EinArray{T}, ind) where {T} = map_prod(A.xs, ind, A.x_indexers)
@inline Base.getindex(A::EinArray{T}, inds::Int...) where {T} = map_prod(A.xs, inds, A.x_indexers)

using StaticArrays:MVector
using Base.Cartesian
using BitBasis

function chasing_game(g, f, xs::NTuple{2}) where N
    xa, xb = xs
    la, lb = 1, 1
    while lb <= length(xb) && la <= length(xa)
        fa = f(xa[la])
        fb = f(xb[lb])
        while fa != fb
            if fa < fb
                la += 1
                la > length(xa) && return
                fa = f(xa[la])
            else
                lb += 1
                lb > length(xb) && return
                fb = f(xb[lb])
            end
        end
        # get number of valid a
        na = 0
        while la+na <= length(xa) && f(xa[la+na]) == fb
            na += 1
        end

        nb = 0
        while lb+nb <= length(xb) && f(xb[lb+nb]) == fa
            nb += 1
        end
        for ka=la:la+na-1, kb=lb:lb+nb-1
            g(ka, kb)
        end
        la += na
        lb += nb
    end
end

function naive_chase(a::Vector{T}, b::Vector{T}) where T
    res = Tuple{Int,Int}[]
    for ia in 1:length(a), ib in 1:length(b)
        a[ia] == b[ib] && push!(res, (ia, ib))
    end
    return res
end

function chasing_game(xs::NTuple)
    res = Tuple{Int,Int}[]
    chasing_game((x,y)->push!(res, (x,y)), x->x, xs)
    return res
end

using Test
using Random
@testset "chasing gate 2 tensors" begin
    cg = chasing_game(([1,2,3], [2,3,4]))
    @test cg == [(2,1), (3,2)]
    cg = chasing_game(([1,2,2,3], [3,4]))
    @test cg == [(4,1)]
    cg = chasing_game(([1,2,2,3], [3,3,4]))
    @test cg == [(4,1), (4,2)]
    cg = chasing_game(([1,2,2,3,3], [3,3,4]))
    @test cg == [(4,1), (4,2), (5,1), (5,2)]
    @test cg == naive_chase([1,2,2,3,3], [3,3,4])
end

function get_inner(::Val{Ni}, x::BitStr{N,T}) where {N, Ni, T}
    BitStr{Ni,T}(x >> (N-Ni))
end

function get_outer(::Val{Ni}, x::BitStr{N,T}) where {Ni,N,T}
    BitStr{N-Ni,T}(x & bmask(1:N-Ni))
end

function get_outer(ni::Val{Ni}, xs::BitStr...) where {Ni}
    bcat(get_outer.(ni, xs)...)
end

@testset "get inner, outer" begin
    @test get_inner(Val(2), bit"0001111") === bit"00"
    @test get_outer(Val(2), bit"0001111") === bit"01111"
    @test get_outer(Val(2), bit"0001111", bit"001") === bit"011111"
end

function sparse_contract(ni::Val{Ni}, a::BinarySparseTensor{T1,Ti,M}, b::BinarySparseTensor{T2,Ti,N}) where {T1,T2,Ni,N,M,Ti}
    out = get_output_array((a,b), (fill(2, M+N-2Ni)...,))
    ia, va = findnz(a)
    ib, vb = findnz(b)
    #chasing_game(x->get_inner(ni, x), (ia,ib)) do la, lb
    for (la, lb) in naive_chase(get_inner.(ni, a), get_inner.(ni, b))
        @show la, lb
        inda, vala = ia[la], va[la]
        indb, valb = ib[lb], vb[lb]
        indout = get_outer(ni, inda, indb)
        out[indout] += vala*valb
    end
    return out
end

@testset "sparse contract" begin
    Random.seed!(2)
    ta = bstrand(4, 0.2)
    tb = bstrand(4, 0.2)
    TA, TB = Array(ta), Array(tb)
    @show sparse_contract(Val(2), ta, tb)
    @test sum(sparse_contract(Val(2), ta, tb)) ≈ sum(ein"ijkl,ijmn->klmn"(TA, TB))
    @test Array(sparse_contract(Val(2), ta, tb)) ≈ ein"ijkl,ijmn->klmn"(TA, TB)
end

# define einsum for both PairWise and PTrace with BinarySparseTensor to have those operations
# dispatch to loop_einsum, since the default dispatch does not support BinarySparseTensor yet
function einsum(::PairWise, code::EinCode{ixs, iy},
            xs::NTuple{NT,BinarySparseTensor{<:OMEinsum.BLASTPYES}},
            size_dict) where {ixs, iy, NT}
    loop_einsum(code, xs, size_dict)
end

function einsum(::PTrace, code::EinCode{ixs, iy},
            xs::NTuple{NT,BinarySparseTensor{<:OMEinsum.BLASTPYES}},
            size_dict) where {ixs, iy, NT}
    loop_einsum(code, xs, size_dict)
end

@eval function OMEinsum.einsum(::OMEinsum.BatchedContract, code::EinCode, xs::NTuple{NT, BinarySparseTensor{<:OMEinsum.BLASTPYES}}, size_dict) where {NT}
    loop_einsum(code, xs, size_dict)
end

@testset "einsum" begin
    sv = SparseVector([1.0,0,0,1,1,0,0,0])
    t1 = bst(sv)
    t2 = bst(sv)
    T1 = Array(t1)
    T2 = Array(t2)
    @test ein"ijk,jkl->il"(t1,t2) ≈ ein"ijk,jkl->il"(T1,T2)
end

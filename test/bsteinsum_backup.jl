using OMEinsum
include("SparseTensor.jl")

OMEinsum.asarray(x::Number, arr::BinarySparseTensor) where T = BinarySparseTensor(asarray(x))

function get_output_array(xs::NTuple{N, BinarySparseTensor}, size) where N
    out = bst_zeros(promote_type(map(eltype,xs)...), size)
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
struct ChasingGame{N,Ti,FT}
    f::FT
    xs::NTuple{N,Vector{Ti}}
end
ChasingGame(xs::NTuple) = ChasingGame(x->x, xs)

function Base.iterate(cg::ChasingGame{2,Ti}, state=(Ti(1), Ti(1))) where Ti
    la, lb = state
    xa, xb = cg.xs
    while lb <= length(xb)
        fb = cg.f(xb[lb])
        # move a
        while la <= length(xa) && cg.f(xa[la]) < fb
            la += 1
        end
        if la <= length(xa)
            return cg.f(xa[la]) == fb && return (la, lb), (la, lb+1)
        else
            return
        end

        fa = cg.f(xa[la])
        while lb <= length(xb) && cg.f(xb[lb]) < fa
            lb += 1
        end
        lb <= length(xb) && cg.f(xb[lb]) == fa && return (la, lb), (la+1, lb)
    end
end

using Test
@testset "chasing gate 2 tensors" begin
    cg = ChasingGame(([1,2,3], [2,3,4]))
    @test [cg...] == [(2,1), (3,2)]
    cg = ChasingGame(([1,2,2,3], [3,4]))
    @test [cg...] == [(4,1)]
    cg = ChasingGame(([1,2,2,3], [3,3,4]))
    @test [cg...] == [(4,1), (4,2)]
    cg = ChasingGame(([1,2,2,3,3], [3,3,4]))
    @test [cg...] == [(4,1), (4,2), (5,1), (5,2)]
end

#=
function Base.iterate(cg::ChasingGame{N,Ti}, state::Tuple{MVector{N,Int},Ti}=(MVector(ones(Int,N)...),Ti(0))) where {N,Ti}
    locs, frontier = state
    while true
        x = cg.xs[i]
        iloc = locs[i]
        iloc, frontier, signal = find_frontier(x->x, x, iloc, frontier)
        locs[i] = iloc
        if signal
            loc[i] += 1
        else
            i += 1
        end
        signal && return (locs, frontier), (locs, frontier)
        if i==N
            # return valid
    end
end
=#

Base.eltype(cg::ChasingGame{N,Ti}) where {N,Ti} = Tuple{MVector{N,Int},Ti}
Base.eltype(::Type{ChasingGame{N,Ti}}) where {N,Ti} = Tuple{MVector{N,Int},Ti}

# if find front, return (loc of old frontier+1, old frontier, true)
# if does not find front, return (loc of new frontier, new frontier, false)
# if does not find front and overflow, return (a number > N, old frontier, false)
@inline function find_frontier(f, x, loc, front)
    val = f(x[loc])
    if val == front
        return loc, front, true
    elseif val > front
        return loc, val, false
    elseif loc<length(x)  # can't find anything
        return find_frontier(f, x, loc+1, front)
    else
        return loc+1, front, false
    end
end

function naive_chase(a::Vector{T}, b::Vector{T}) where T
    res = T[]
    for ia in a
        if ia in b
            push!(res, ia)
        end
    end
    return res
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

function sparse_contract(a::BinarySparseTensor{T,M}, b::BinarySparseTensor{T,N}, ni::Val{Ni}) where Ni
    cg = ChasingGame(x->get_inner(x, ni), (a,b))
    for ((i,j), match) in cg
        oind = get_outer(ni, (a,b))
    end
end

using Test

@testset "chasing" begin
    @test find_frontier(x->x, [1,2,3], 1, 3) == (4,3,true)
    @test find_frontier(x->x, [1,2,3], 1, 9) == (4,9,false)
    @test find_frontier(x->x, [1,2,7], 1, 3) == (3,7,false)
    a = rand(1:1000, 100) |> sort
    b = rand(1:1000, 100) |> sort
    cg = ChasingGame((a, b))
    res = naive_chase(a,b)
    @test getindex.([cg...], 2) == unique(res)
    a = rand(1:10, 100) |> sort
    b = rand(1:10, 100) |> sort
    cg = ChasingGame((a, b))
    res = naive_chase(a,b)
    @test getindex.([cg...], 2) == unique(res)
end

for x in ChasingGame
    println(x)
end

function chasing_game(g, f, xs::NTuple{2}) where N
    xa, xb = xs
    la, lb = 1, 1
    while lb <= length(xb)
        @inbounds fb = f(xb[lb])
        @inbounds while la <= length(xa) && f(xa[la]) < fb
            la += 1
        end
        la > length(xa) && return
        @inbounds fa = f(xa[la])
        for k=lb:length(xb)
            @inbounds if f(xb[k]) == fa
                g(la, k)
            else
                break
            end
        end
        la += 1
    end
end

function get_inner_indices(A::EinArray, y)
    t1, t2 = A.xs
    Ni = ndims(A.ICIS)
    d = Dict{Ti,Tv}()
    for (i,t) in enumerate(A.ts)
        nzind, nzval = findnz(t)
        iind = readbit(nzind, innerloc...)
        if i==0
        d[iind] = get(d, iind, 1)*nzval
    end
    for iind in inner_indices(A)
        # index propagation
        res = matchindices(ind, A.xs[2:end])
        for ri in res
            y[ri.yind] = prod(ri.xvals)
        end
    end
end

function Base._mapreducedim!(f, op, y::BinarySparseTensor, A::EinArray)
    A
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

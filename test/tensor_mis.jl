using MISAlgorithms
using OMEinsum, Test

g1 = EliminateGraph(5, [1=>2, 2=>3, 2=>4, 3=>4, 4=>5])

@test mis2(g1) == 3

struct InfNum <: Number
    n::Rational{Int}
end

∅ = InfNum(-1)

const e = InfNum(0)
const ∞ = InfNum(1)

function Base.show(io::IO, inf::InfNum)
    print(io,"∞($(inf == ∅ ? "∅" : inf.n))")
end

#Base.:^(a::InfNum, b::Int) = a == ∅ ? a : InfNum(a.n*b)
Base.:*(a::InfNum, b::InfNum) = a == ∅ || b == ∅ ? ∅ : InfNum(a.n+b.n)
Base.:+(a::InfNum, b::InfNum) = InfNum(max(a.n,b.n))
Base.zero(::Type{InfNum}) = ∅

T(a::Rational, b::Rational) = [e InfNum(b); InfNum(a) ∅]

# show the final result
res = ein"ab,bc,cd,bd,de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
res = ein"(((ab,bc),cd),bd),de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
@test res == 3

# show the full table
tbl = ein"ab,bc,cd,bd,de->abcde"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))

isinferier(tbl, σA::NTuple{X,Int}, σB::NTuple{X,Int}) where X = σA!=σB && tbl[σA...].n <= tbl[σB...].n && all(σA .>= σB)
isinferier(tbl, σA::NTuple{X,Int}) where X = σA == ∅ || any(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))

# this marks eliminatable elements
inferier_table(tbl, σA::NTuple{X,Int}) where X = map(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))
inferier_table(tbl) where X = map(σB -> isinferier(tbl, σB.I), CartesianIndices(tbl))

@test !isinferier(tbl, (1,1,1,1,1), (1,1,1,1,1))
@test isinferier(tbl, (1,2,2,2,1), (1,1,1,1,1))
@test isinferier(tbl, (1,2,2,2,1))
@test !isinferier(tbl, (1,1,1,1,1))

inferier_table(tbl, (1,1,1,1,1))
inferier_table(tbl)

@test sum(tbl .!= ∅) == 12
@test sum(inferier_table(tbl)) == 32 - 12

# show the part table, (a, b) are outer legs
tbl = ein"ab,bc,cd,bd,de->ab"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))

@test sum(tbl .!= ∅) == 3
@test isinferier(tbl, (1,2), (1,1))
@test sum(inferier_table(tbl)) == 2

# sub-graph: triangles
# 1.2.1
ein"ij,jk,ki->jk"(T(1//2,1//2), T(1//2,1//2), T(1//2,1//2))

tsplib, satlib, qbflib, junction tree
# 1.2.2.1
ein"ij,ik,jl,kl->l"(T(1//2,1//2), T(1//2,1//2), T(1//2,1//2), T(1//2,1//2))
# 1.2.2.2
ein"ij,ik->jk"(T(1//2,1//1), T(1//2,1//1))
ein"ij,ik,jl,kl,km->lm"(T(1//2,1//2), T(1//2,1//3), T(1//2,1//2), T(1//3,1//2), T(1//3,1//1))
ein"ij,ik,jl,kn,km->lm"(T(1//2,1//2), T(1//2,1//3), T(1//2,1//1), T(1//3,1//1), T(1//3,1//1))
ein"ij,ik,jl,km->lm"(T(1//2,1//2), T(1//2,1//2), T(1//2,1//1), T(1//2,1//1))

# graph K33 + 1, verify the mirror rule: include v, or eliminate v and its mirrors.
T33 = T(1//3, 1//3)
T34 = T(1//3, 1//4)
ein"ia,ib,ic,ja,jb,jc,ka,kb,kc,ab->ijk"(T34, T34, T33,T34,T34,T33,T34,T34,T33,T(1//4,1//4))

# neighbor analysis
# degree == 2, from 1.2.2.2 and 1.2.1, we know the maximum branching is 2, <= 2^(1/3).
#

# TODO:
# 1. sparse tensor type to store possible configurations
# 2. merging two tensors, a flaw of this context un-aware algebra?
#    maybe it is also an advantage, it represents probability, gives a greedy version of algorithm.

# not working now, since gcd is not properly implemented for GPU, should be an easy fix.
#using CuArrays
#ein"ij,jk,ki->jk"(CuArray(T(1//2,1//2)), CuArray(T(1//2,1//2)), CuArray(T(1//2,1//2)))

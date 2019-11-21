using MISAlgorithms
using OMEinsum, Test

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

isinferier(tbl, σA::NTuple{X,Int}, σB::NTuple{X,Int}) where X = σA!=σB && tbl[σA...].n <= tbl[σB...].n && all(σA .>= σB)
isinferier(tbl, σA::NTuple{X,Int}) where X = σA == ∅ || any(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))

# this marks eliminatable elements
inferier_table(tbl, σA::NTuple{X,Int}) where X = map(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))
inferier_table(tbl) where X = map(σB -> isinferier(tbl, σB.I), CartesianIndices(tbl))

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

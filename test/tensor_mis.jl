using MISAlgorithms
using OMEinsum, Test

g1 = EliminateGraph(5, [1=>2, 2=>3, 2=>4, 3=>4, 4=>5])

@test mis2(g1) == 3

struct InfNum <: Number
    n::Rational
end

∅ = InfNum(-1)

const e = InfNum(0)
const ∞ = InfNum(1)

function Base.show(io::IO, inf::InfNum)
    print(io, inf == ∅ ? "∅" : "∞($(inf.n))")
end

#Base.:^(a::InfNum, b::Int) = a == ∅ ? a : InfNum(a.n*b)
Base.:*(a::InfNum, b::InfNum) = a == ∅ || b == ∅ ? ∅ : InfNum(a.n+b.n)
Base.:+(a::InfNum, b::InfNum) = InfNum(max(a.n,b.n))
Base.zero(::Type{InfNum}) = ∅

T(a::Rational, b::Rational) = [e InfNum(b); InfNum(a) ∅]

# show the final result
res = ein"ab,bc,cd,bd,de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
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

isinferier(tbl, (1,2,2,2,1), (1,1,1,1,1))
inferier_table(tbl, (1,1,1,1,1))
inferier_table(tbl)

@test sum(tbl .!= ∅) == 12
@test sum(inferier_table(tbl)) == 32 - 12

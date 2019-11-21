include("tensor_mis.jl")

ST(a::Rational, b::Rational) = BinarySparseTensor(SparseVector{InfNum,Int}(vec([e InfNum(b); InfNum(a) ∅])))
Base.iszero(a::InfNum) = a==∅
Base.:(==)(a::InfNum, b::Int) = b==0 ? iszero(a) : throw(MethodError(:(==), (a,b)))

A = ST(1//2, 1//3)
@test nnz(A) == 3
SparseArrays.dropzeros(A::BinarySparseTensor) = BinarySparseTensor(dropzeros(A.data))
dropzeros(A)

function drop_bad(A::BinarySparseTensor{InfNum})
    mask = ones(Bool, nnz(A))
    for (i,(si, sv)) = enumerate(findnz(A))
        if isinferier(A, si)
            mask[i] = false
        end
    end
    sv = SparseVector(sv)
    A = BinarySparseTensor(sv)
end

# show the final result
res = ein"ab,bc,cd,bd,de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
res = ein"(((ab,bc),cd),bd),de->"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))[].n
@test res == 3

# show the full table
tbl = ein"ab,bc,cd,bd,de->abcde"(T(1//1,1//3),T(1//3,1//2),T(1//2,1//3),T(1//3,1//3),T(1//3,1//1))

"""
1 if a bit in a is larger or equal to the bit b.
"""
⪰(a::T, b::T) where {N,T<:BitStr{N}} = (a | ~b) & bmask(1:N)
≻(a::T, b::T) where {N,T<:BitStr{N}} = (a & ~b) & bmask(1:N)

(bit"1010" ⪰ bit"0110") == bit"1011"
(bit"1010" ≻ bit"0110") == bit"1000"

"""
1 if a bit in a is smaller or equal to the bit b.
"""
⪯(a::T, b::T) where {N,T<:BitStr{N}} = (~a | b) & bmask(1:N)
≺(a::T, b::T) where {N,T<:BitStr{N}} = (~a & b) & bmask(1:N)

(bit"1010" ⪯ bit"0110") == bit"0111"
(bit"1010" ≺ bit"0110") == bit"0100"

isinferier(tbl::BinarySparseTensor, σA::BitStr, σB::BitStr) where X = σA!=σB && tbl[σA...].n <= tbl[σB...].n && all(σA .>= σB)
isinferier(tbl, σA::NTuple{X,Int}) where X = σA == ∅ || any(σB -> isinferier(tbl, σA, σB.I), CartesianIndices(tbl))

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

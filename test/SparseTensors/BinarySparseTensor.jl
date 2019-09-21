using BitBasis, SparseArrays, MISAlgorithms
using Test, Random

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
    AT = Array(t)
    @test permutedims(t, (2,1,3)) == permutedims(AT, (2,1,3))
    @test t == AT

    m = zeros(Int, 2,2,2); m[[1,4,5]] .= 1
    @test Array(t) == m
    @test vec(t) == vec(m)
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
    @test MISAlgorithms.bpermute(0b1100, [2,4,3,1]) === 0b0110
    @test MISAlgorithms.bpermute(0b1100, [4,2,3,1]) === 0b0101
    @test bpermute(bit"1100", [2,4,3,1]) === bit"0110"
    @test bpermute(bit"1100", [4,2,3,1]) === bit"0101"
    @test sort(BitStr64{5}[1,8,4,2,9]) == sort([1,8,4,2,9])
end

@testset "bst einsum" begin
    include("bsteinsum.jl")
end

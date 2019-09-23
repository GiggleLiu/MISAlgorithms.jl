using OMEinsum, MISAlgorithms, BitBasis
using MISAlgorithms: get_outer, get_inner, get_inner_and_batch, get_batch, chasing_game
using Test
using Random, SparseArrays

function naive_chase(a::Vector{T}, b::Vector{T}) where T
    res = Tuple{Int,Int}[]
    for ia in 1:length(a), ib in 1:length(b)
        a[ia] == b[ib] && push!(res, (ia, ib))
    end
    return res
end

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

@testset "get inner, outer" begin
    # indices are "iiiooobbb"
    @test get_inner(Val(2), bit"0001111") === bit"00"
    @test get_inner_and_batch(Val(2), Val(1), bit"0001111") === bit"000"
    @test get_batch(Val(2), Val(1), bit"0001111") === bit"0"
    @test get_outer(Val(2), Val(0), bit"0001111") === bit"01111"
    @test get_outer(Val(2), Val(0), bit"0001111", bit"001") === bit"101111"
    @test get_outer(Val(2), Val(1), bit"0001111", bit"000") === bit"01111"
end

@testset "sparse contract" begin
    Random.seed!(2)
    ta = bstrand(4, 0.5)
    tb = bstrand(4, 0.5)
    TA, TB = Array(ta), Array(tb)
    @test sum(sparse_contract(Val(2), Val(0), ta, tb)) ≈ sum(ein"lkji,nmji->lknm"(TA, TB))
    @test Array(sparse_contract(Val(2), Val(0), ta, tb)) ≈ ein"lkji,nmji->lknm"(TA, TB)

    # batched
    ta = bstrand(5, 0.5)
    tb = bstrand(5, 0.5)
    TA, TB = Array(ta), Array(tb)
    @test sum(Array(sparse_contract(Val(2), Val(1), ta, tb))) ≈ sum(ein"lkbji,nmbji->lknmb"(TA, TB))
    @test Array(sparse_contract(Val(2), Val(1), ta, tb)) ≈ ein"lkbji,nmbji->lknmb"(TA, TB)
end

@testset "einsum batched contract" begin
    Random.seed!(2)

    perms = MISAlgorithms.analyse_batched_perm(('b','j','c','a','e'), ('k','d','c','a','e'), ('b','j','k','d','c'))
    @test perms == ([1, 2, 3, 4, 5], [1, 2, 3, 4, 5], (1, 2, 3, 4, 5), 2, 1)

    sv = SparseVector([1.0,0,0,1,1,0,0,0])
    t1 = bst(sv)
    t2 = bst(sv)
    T1 = Array(t1)
    T2 = Array(t2)
    @test ein"ijk,jkl->il"(t1,t2) ≈ ein"ijk,jkl->il"(T1,T2)

    ta = bstrand(2, 0.5)
    tb = bstrand(2, 0.5)
    TA, TB = Array(ta), Array(tb)
    @test ein"ij,jk->ik"(ta,tb) ≈ ein"ij,jk->ik"(TA,TB)
    @test ta == TA
    @test tb == TB

    ta = bstrand(7, 0.5)
    tb = bstrand(6, 0.5)
    TA, TB = Array(ta), Array(tb)
    code = ein"ijklmbc,ijbxcy->bcmlxky"
    @test OMEinsum.match_rule(code) == OMEinsum.BatchedContract()
    res = code(ta, tb)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ code(TA, TB)
end

@testset "sum, ptrace and permute" begin
    ta = bstrand(7, 0.7)
    TA = Array(ta)
    res = dropsum(ta, dims=(2,4))
    @test res isa BinarySparseTensor
    @test Array(res) ≈ dropsum(TA, dims=(2,4))
    @test dropsum(ta) ≈ dropsum(TA)
    @test sum(ta) ≈ sum(TA)
    # sum
    res = ein"ijklbca->"(ta)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ ein"ijklbca->"(TA)
    res = ein"ijklbca->i"(ta)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ ein"ijklbca->i"(TA)

    # trace
    tb = bstrand(2, 1.0)
    TB = Array(tb)
    res = ein"ii->"(tb)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ ein"ii->"(TB)
    @test OMEinsum.match_rule(ein"ijjlbca->ailcb") == OMEinsum.PTrace()

    # reduction
    code = ein"ijjlbbb->ijlb"
    res = code(ta)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ code(TA)
    code = ein"ijjlbbb->ljbi"
    res = code(ta)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ code(TA)

    # ptrace
    res = ein"ijjlbca->ailcb"(ta)
    @test res isa BinarySparseTensor
    @test Array(res) |> sum ≈ ein"ijjlbca->ailcb"(TA) |> sum
    @test Array(res) ≈ ein"ijjlbca->ailcb"(TA)

    # permute
    res = ein"ijklbca->abcijkl"(ta)
    @test res isa BinarySparseTensor
    @test Array(res) ≈ ein"ijklbca->abcijkl"(TA)
end

@testset "clean up tensors" begin
    ta = bstrand(7, 0.5)
    TA = Array(ta)
    # first reduce indices
    for code in [ein"iiiiiii->iiiiii", ein"iikjjjj->ikj", ein"iikjjjl->ikj"]
        res = code(ta)
        @test res isa BinarySparseTensor
        @test res ≈ code(TA)
    end
end

@testset "clean up" begin
    @test MISAlgorithms._ymask_from_reds(Int, 5, [[2,3], [4,1]]) == bmask(5,2,4)
    @test !OMEinsum.match_rule(MISAlgorithms.IndexReduction(), ein"ijk->ijk")
    @test OMEinsum.match_rule(MISAlgorithms.IndexReduction(), ein"ijkj->ijk")
    @test !OMEinsum.match_rule(MISAlgorithms.IndexReduction(), ein"ijkjl->ijk")

    @test allsame(bit"000110", bmask(BitStr64{6}, 2,3))
    @test !allsame(bit"000110", bmask(BitStr64{6}, 2,4))
    @test MISAlgorithms._get_reductions((1,2,2,4,3,1,5), (1,2,3,4,5)) == [[1,6], [2,3]]

    @test MISAlgorithms.getalllabels(ein"ijk,jkl->oo") == ['i', 'j', 'k', 'l', 'o']
end


@testset "count legs" begin
    @test MISAlgorithms.count_legs((1,2), (2,3), (1,3)) == Dict(1=>2,2=>2,3=>2)
    @test MISAlgorithms.dangling_legs(((1,1,2), (2,3)), (5,7)) == (((), (2,)), (1,2))
    @test MISAlgorithms.dangling_nleg_labels(((1,1,2), (2,3)), (5,7)) == (((1,), (3,)), (5,7))
end

MISAlgorithms.unset(bit"111100", bit"001101") == bit"110000"

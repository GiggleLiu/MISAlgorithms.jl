OMEinsum.asarray(x::Number, arr::BinarySparseTensor) where T = BinarySparseTensor(asarray(x))

function OMEinsum.get_output_array(xs::NTuple{N, BinarySparseTensor{Tv,Ti,M} where {Tv,M}}, size) where {N,Ti}
    out = bst_zeros(promote_type(map(eltype,xs)...), Ti, length(size))
end

function chasing_game(g, f, xs::NTuple{2,Vector}) where N
    xa, xb = xs
    la, lb = 1, 1
    while lb <= length(xb) && la <= length(xa)
        @inbounds fa = f(xa[la])
        @inbounds fb = f(xb[lb])
        @inbounds while fa != fb
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
        @inbounds while la+na <= length(xa) && f(xa[la+na]) == fb
            na += 1
        end

        nb = 0
        @inbounds while lb+nb <= length(xb) && f(xb[lb+nb]) == fa
            nb += 1
        end
        for ka=la:la+na-1, kb=lb:lb+nb-1
            g(ka, kb)
        end
        la += na
        lb += nb
    end
end

function chasing_game(xs::NTuple)
    res = Tuple{Int,Int}[]
    chasing_game((x,y)->push!(res, (x,y)), x->x, xs)
    return res
end

@generated function ibcat(bits::NTuple{N, BitStr{M,T} where M}) where {N,T}
    total_bits = BitBasis.sum_length(bits.parameters...)

    quote
        val, len = T(0), 0
        @nexprs $N k->(val += buffer(bits[k]) << len; len += length(bits[k]))
        return BitStr{$total_bits,T}(val)
    end
end

function get_inner(::Val{Ni}, x::BitStr{N,T}) where {N, Ni, T}
    BitStr{Ni,T}(x >> (N-Ni))
end

function get_batch(::Val{Ni}, ::Val{Nb}, x::BitStr{N,T}) where {N, Ni, Nb, T}
    BitStr{Nb,T}((x >> (N-Ni-Nb)) & bmask(1:Nb))
end

function get_inner_and_batch(::Val{Ni}, ::Val{Nb}, x::BitStr{N,T}) where {N, Nb, Ni, T}
    BitStr{Nb+Ni,T}(x >> (N-Nb-Ni))
end

function get_outer(::Val{Ni}, ::Val{Nb}, x::BitStr{N,T}) where {Ni,Nb,N,T}
    BitStr{N-Ni-Nb,T}(x & bmask(1:N-Ni-Nb))
end

function get_outer(ni::Val{Ni}, nb::Val{Nb}, xs::BitStr...) where {Ni,Nb}
    ibcat((get_outer.(ni, nb, xs)...,get_batch(ni, nb, xs[1])))
end

function sparse_contract(ni::Val{Ni}, nb::Val{Nb}, a::BinarySparseTensor{T1,Ti,M}, b::BinarySparseTensor{T2,Ti,N}) where {T1,T2,Ni,Nb,N,M,Ti}
    out = OMEinsum.get_output_array((a,b), (fill(2, M+N-2Ni-Nb)...,))
    ia, va = findnz(a)
    ib, vb = findnz(b)
    chasing_game(x->get_inner_and_batch(ni, nb, x), (ia,ib)) do la, lb
    # for (la, lb) in naive_chase(get_inner_and_batch.(ni, nb, ia), get_inner.(ni, nb, ib))
        inda, vala = ia[la], va[la]
        indb, valb = ib[lb], vb[lb]
        indout = get_outer(ni, nb, inda, indb)
        out[indout] += vala*valb
    end
    return out
end

# can be used in either static or dynamic invoke
function analyse_batched_perm(iAs, iBs, iOuts)
    iABs = iAs ∩ iBs
    pres   = iABs ∩ iOuts
    broad  = setdiff((iAs ∩ iOuts) ∪ (iBs ∩ iOuts), pres)
    summed = setdiff(iABs, pres)

    iAps, iAbs, iAss = pres ∩ iAs, broad ∩ iAs, summed ∩ iAs
    iBps, iBbs, iBss = pres ∩ iBs, broad ∩ iBs, summed ∩ iBs

    pA   = OMEinsum.indexpos.(Ref(iAs), vcat(iAbs, iAps, iAss))
    pB   = OMEinsum.indexpos.(Ref(iBs), vcat(iBbs, iBps, iBss))
    iABs = vcat(iAbs, iBbs, iAps)
    pOut = OMEinsum.indexpos.(Ref(iABs), iOuts)

    return pA, pB, pOut, length(iAss), length(iAps)
end

@eval function OMEinsum.einsum(::OMEinsum.BatchedContract, code::EinCode{ixs, iy}, xs::NTuple{NT, BinarySparseTensor}, size_dict) where {NT, ixs, iy}
    a, b = xs
    pa, pb, pout, Ni, Nb = analyse_batched_perm(ixs..., iy)
    A = copy(a)
    B = copy(b)
    a = permutedims(a, pa)
    b = permutedims(b, pb)
    out = sparse_contract(Val(Ni), Val(Nb), a, b)
    permutedims(out, pout)
end

function OMEinsum.einsum(sm::OMEinsum.EinRule, code::EinCode{ixs, iy}, xs::NTuple{NT, BinarySparseTensor}, size_dict) where {ixs, iy, NT}
    throw(ArgumentError("Eincode $code not supported for BinarySparseTensor yet."))
end
function OMEinsum.einsum(sm::OMEinsum.MatMul, code::EinCode{ixs, iy}, xs::NTuple{NT, BinarySparseTensor}, size_dict) where {ixs, iy, NT}
    einsum(OMEinsum.BatchedContract(), code, xs, size_dict)
end

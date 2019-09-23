export allin
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

############ BitBasis
@generated function ibcat(bits::NTuple{N, BitStr{M,T} where M}) where {N,T}
    total_bits = BitBasis.sum_length(bits.parameters...)

    quote
        val, len = T(0), 0
        @nexprs $N k->(val += buffer(bits[k]) << len; len += length(bits[k]))
        return BitStr{$total_bits,T}(val)
    end
end

export allsame
allsame(x::T, mask::T) where T<:Integer = allone(x, mask) || !anyone(x, mask)

######################### index manipulation
function count_legs(ixs::(NTuple{N,T} where N)...) where {T}
    lc = Dict{T,Int}()
    for l in Iterators.flatten(ixs)
        lc[l] = get(lc, l, 0) + 1
    end
    return lc
end

"""
return positions of dangling legs.
"""
function dangling_legs(ixs, iy, lc=count_legs(ixs..., iy))
    _dlegs.(ixs, Ref(lc)), _dlegs(iy, lc)
end
_dlegs(ix, lc) = (findall(iix->lc[iix] == 1, [ix...])...,)

"""
return positions of dangling legs.
"""
function dangling_nleg_labels(ixs, iy, lc=count_legs(ixs..., iy))
    _dnlegs.(ixs, Ref(lc)), _dnlegs(iy, lc)
end
function _dnlegs(ix, lc)
    (unique(filter(iix->count(==(iix), ix)==lc[iix], [ix...]))...,)
end

function dumplicated_legs(ix)
    labels = OMEinsum.tunique(ix)
    findall(l->count(==(l), ix)>1, labels)
end

function getalllabels(code::EinCode{ixs, iy}) where {ixs, iy}
    reduce(∪, (ixs..., iy))
end

function allin(x, y)
    all(xi->xi in y, x)
end

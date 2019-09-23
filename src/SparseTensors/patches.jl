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

@generated function ibcat(bits::NTuple{N, BitStr{M,T} where M}) where {N,T}
    total_bits = BitBasis.sum_length(bits.parameters...)

    quote
        val, len = T(0), 0
        @nexprs $N k->(val += buffer(bits[k]) << len; len += length(bits[k]))
        return BitStr{$total_bits,T}(val)
    end
end

######################### index manipulation
function count_legs(ixs::(NTuple{N,T} where N)...) where {T}
    keys, vals = T[], Int[]
    for l in Iterators.flatten(ixs)
        res = findfirst(==(l), keys)
        if res === nothing
            push!(keys, l)
            push!(vals, 1)
        else
            vals[res] += 1
        end
    end
    return LEGCOUNT{(keys...,), (vals...,)}
end

struct LEGCOUNT{lbs, ns} end
function Base.getindex(::Type{LEGCOUNT{lbs, ns}}, lb) where {lbs, ns}
    ns[findfirst(==(lb), lbs)]
end

"""
return positions of dangling legs.
"""
function dangling_legs(ixs, iy, lc::Type{<:LEGCOUNT}=count_legs(ixs, iy))
    _dlegs.(ixs, lc), _dlegs(iy, lc)
end
_dlegs(ix, lc) = (findall(iix->lc[iix] == 1, [ix...])...,)

function dumplicated_legs(ix)
    labels = OMEinsum.tunique(ix)
    findall(l->count(==(l), ix)>1, labels)
end

using Test
lc = MISAlgorithms.LEGCOUNT{(1,2,3), (6,7,8)}
@test lc[2] == 7
@test MISAlgorithms.count_legs(((1,2), (2,3)), (1,3)) == LEGCOUNT{(1,2,3), (2,2,2)}
@test MISAlgorithms.dangling_legs(((1,2), (2,3)), (5,7)) == (((1,), (2,)), (1,2))

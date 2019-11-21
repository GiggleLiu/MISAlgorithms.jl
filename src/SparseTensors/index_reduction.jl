export dropsum

using TupleTools
using OMEinsum: allunique, tunique

function OMEinsum.einsum(sm::OMEinsum.Sum, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    dims = (findall(i -> i ∉ iy, ixs[1])...,)
    (ix1,) = ixs
    ix1f = filter!(i -> i in iy, collect(ix1))
    perm = map(i -> findfirst(==(i), ix1f), iy)
    res = dropsum(xs[1], dims=dims)
    perm == iy ? res : permutedims(res, perm)
end

function cleanup_dangling_nlegs(ixs, xs, iy)
    lc = count_legs(ixs..., iy)
    danglegsin, danglegsout = dangling_nleg_labels(ixs, iy, lc)
    newxs = Any[xs...]
    newixs = Any[ixs...]
    newiy = iy
    for i = 1:length(xs)
        ix, dlx = newixs[i], danglegsin[i]
        if !isempty(dlx)
            newxs[i] = multidropsum(xs[i], dims=[findall(==(l), ix) for l in dlx])
            newixs[i] = (filter(l->l∉dlx, [ix...])...,)
        end
    end
    if !isempty(danglegsout)
        newiy = ([l for l in newiy if l ∉ danglegsout]...,)
    end
    return newixs, newxs, newiy
end

function cleanup_dumplicated_legs(ixs, xs, iy)
    newxs = Any[xs...]
    newixs = Any[ixs...]
    for i in 1:length(xs)
        ix = ixs[i]
        if !allunique(ix)
            newix = (tunique(ix)...,)
            newxs[i] = einsum(IndexReduction(), EinCode{(ix,), newix}(), (xs[i],), OMEinsum.get_size_dict((ix,), (xs[i],)))
            newixs[i] = newix
        end
    end
    newiy = iy
    if !allunique(newiy)
        newiy = (tunique(newiy)...,)
    end
    return newixs, newxs, newiy
end

function OMEinsum.einsum(sm::OMEinsum.DefaultRule, code::EinCode{ixs, iy}, xs::NTuple{NT, BinarySparseTensor}, size_dict) where {ixs, iy, NT}
    if OMEinsum.match_rule(IndexCopy(), code)
        return einsum(IndexCopy(), code, xs, size_dict)
    elseif OMEinsum.match_rule(IndexBroadcast(), code)
        return einsum(IndexBroadcast(), code, xs, size_dict)
    end
    # dangling legs or multi-legs
    ycode = []
    newixs, newxs, newiy = cleanup_dangling_nlegs(ixs, xs, iy)
    newiy!=iy && pushfirst!(ycode, EinCode{(newiy,), iy}())
    newixs, newxs, newnewiy = cleanup_dumplicated_legs(newixs, newxs, newiy)
    newiy!=newnewiy && pushfirst!(ycode, EinCode{(newnewiy,), newiy}())
    @show ycode
    newixs == ixs && newnewiy == iy && throw(ArgumentError("Einsum not implemented for $code"))
    # dangling merge multi-legs to one leg
    res = einsum(EinCode{(newixs...,), newnewiy}(), newxs, size_dict)
    for code in ycode
        # TODO: broadcast and duplicate
        res = einsum(code, (res,), size_dict)
    end
    return res
end

struct IndexReduction<:OMEinsum.EinRule{1} end
struct IndexCopy<:OMEinsum.EinRule{1} end
struct IndexBroadcast<:OMEinsum.EinRule{1} end

function OMEinsum.match_rule(::IndexReduction, code::EinCode{ixs, iy}) where {ixs, iy}
    length(ixs) > 1 && return false
    ix = ixs[1]
    (allunique(ix) || !allunique(iy)) && return false
    allin(iy, ix) && allin(ix, iy) || return false
    return true
end

function OMEinsum.match_rule(::IndexCopy, code::EinCode{ixs, iy}) where {ixs, iy}
    length(ixs) > 1 && return false
    ix = ixs[1]
    (!allunique(ix) || allunique(iy)) && return false
    allin(iy, ix) && allin(ix, iy) || return false
    return true
end

function OMEinsum.match_rule(::IndexBroadcast, code::EinCode{ixs, iy}) where {ixs, iy}
    length(ixs) > 1 && return false
    ix = ixs[1]
    !allunique(ix) && return false
    allin(ix, iy) || return false
    filter(iiy->(iiy in ix), [iy...]) == [ix...] || return false
    return true
end

function _get_reductions(ix, iy)
    reds = []
    for l in iy
        count(==(l), ix) > 1 && push!(reds, findall(==(l), ix))
    end
    return reds
end

function _get_traces(ix, iy)
    reds = []
    for l in tunique(ix)
        l ∉ iy && count(==(l), ix) > 1 && push!(reds, findall(==(l), ix))
    end
    return reds
end

function OMEinsum.einsum(sm::IndexReduction, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    reduce_indices(xs[1], _get_reductions(ixs[1], iy))
end

function OMEinsum.einsum(sm::IndexCopy, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    ix = ixs[1]
    copy_indices(xs[1], map(l->(findall(==(l), iy)...,), ix))
end

function OMEinsum.einsum(sm::IndexBroadcast, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    ix = ixs[1]
    throw(ArgumentError("Not implemented"))
end

function OMEinsum.einsum(sm::OMEinsum.PTrace, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    res = trace_indices(xs[1]; dims=_get_traces(ixs[1], iy))
    newix = (filter(ix->ix ∈ iy, [ixs[1]...])...,)
    return einsum(EinCode{(newix,), iy}(), (res,), size_dict)
end

function _ymask_from_reds(::Type{Ti}, ndim::Int, reds) where Ti
    ymask = flip(Ti(0), bmask(Ti, 1:ndim))
    for red in reds
        ymask = unsetbit(ymask, bmask(Ti,red[2:end]...))
    end
    return ymask
end

function unsetbit(x::T, mask::T) where T<:Integer
    msk = ~T(0) ⊻ mask
    x & msk
end

function _ymask_from_trs(::Type{Ti}, ndim::Int, reds) where Ti
    ymask = flip(Ti(0), bmask(Ti, 1:ndim))
    for red in reds
        ymask = unsetbit(ymask, bmask(Ti, red...))
    end
    return ymask
end

function reduce_indices(t::BinarySparseTensor{Tv,Ti}, reds) where {Tv,Ti}
    inds = Ti[]
    vals = Tv[]
    ymask = _ymask_from_reds(Ti, ndims(t), reds)
    bits = baddrs(ymask)
    red_masks = [bmask(Ti, red...) for red in reds]
    for (ind, val) in zip(t.data.nzind, t.data.nzval)
        b = ind-1
        if all(red->allsame(b, red), red_masks)
            b = readbit(b, bits...)
            push!(inds, b+1)
            push!(vals, val)
        end
    end
    order = sortperm(inds)
    return BinarySparseTensor(SparseVector(1<<length(bits), inds[order], vals[order]))
end

function copy_indices(t::BinarySparseTensor{Tv,Ti}, targets) where {Tv,Ti}
    inds = Ti[]
    vals = Tv[]
    nbits = sum(length, targets)
    for (ind, val) in zip(t.data.nzind, t.data.nzval)
        b = ind-1
        b = copybits(b, targets)
        push!(inds, b+1)
        push!(vals, val)
    end
    order = sortperm(inds)
    @show nbits
    return BinarySparseTensor(SparseVector(1<<nbits, inds[order], vals[order]))
end

function copybits(b::Ti, targets) where Ti
    res = Ti(0)
    for (i,t) in enumerate(targets)
        for it in t
            res |= readbit(b, i)<<(it-1)
        end
    end
    return res
end

function trace_indices(t::BinarySparseTensor{Tv,Ti}; dims) where {Tv,Ti}
    ymask = _ymask_from_trs(Ti, ndims(t), dims)
    bits = baddrs(ymask)
    red_masks = [bmask(Ti, red...) for red in dims]
    NO = length(bits)
    sv = SparseVector(1<<NO, Ti[], Tv[])
    for (ind, val) in zip(t.data.nzind, t.data.nzval)
        b = ind-1
        if all(red->allsame(b, red), red_masks)
            b = readbit(b, bits...)
            sv[b+1] += val
        end
    end
    return BinarySparseTensor(sv)
end

Base._sum(f, t::BinarySparseTensor, ::Colon) = Base._sum(f, t.data, Colon())
function _dropsum(f, t::BinarySparseTensor{Tv,Ti,N}, dims) where {Tv,Ti,N}
    remdims = (setdiff(1:N, dims)...,)
    Tf = typeof(f(Tv(0)))
    d = Dict{Ti,Tf}()
    for (ind, val) in zip(t.data.nzind, t.data.nzval)
        rd = readbit(ind-1, remdims...)
        d[rd] = get(d, rd, Tf(0)) + f(val)
    end
    ks = collect(keys(d))
    order = sortperm(ks)
    vals = collect(values(d))[order]
    return bst(SparseVector(1<<length(remdims), ks[order].+1, vals))
end

_dropsum(f, t::BinarySparseTensor, dims::Colon) = Base._sum(f, t, dims)
_dropsum(f, t::AbstractArray, dims::Colon) = Base._sum(f, t, dims)
_dropsum(f, t::AbstractArray, dims) = dropdims(Base._sum(f, t, dims), dims=dims)
dropsum(t::AbstractArray; dims=:) = _dropsum(identity, t, dims)

function multidropsum(t::BinarySparseTensor; dims)
    all(d->length(d) == 1, dims) && return dropsum(t; dims=first.(dims))
    trace_indices(t; dims=dims)
end

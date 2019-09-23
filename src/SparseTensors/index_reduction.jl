export dropsum

using TupleTools
function OMEinsum.einsum(sm::OMEinsum.Sum, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    dims = (findall(i -> i ∉ iy, ixs[1])...,)
    (ix1,) = ixs
    ix1f = filter!(i -> i in iy, collect(ix1))
    perm = map(i -> findfirst(==(i), ix1f), iy)
    res = dropsum(xs[1], dims=dims)
    perm == iy ? res : permutedims(res, perm)
end

using OMEinsum: allunique, tunique
function OMEinsum.einsum(sm::OMEinsum.DefaultRule, code::EinCode{ixs, iy}, xs::NTuple{NT, BinarySparseTensor}, size_dict) where {ixs, iy, NT}
    lc = count_legs(ixs..., iy)
    danglegsin, dangdlegsout = dangling_legs(ixs, iy, lc)
    newxs = Any[xs...]
    newixs = Any[ixs...]
    ycode = Any[]
    for i = 1:NT
        ix, dlx = ixs[i], danglegsin[i]
        if !isempty(dlx)
            newxs[i] = dropsum(xs[i], dims=dlx)
            newixs[i] = TupleTools.deleteat(ixs[i], dlx)
        end
    end
    if !isempty(danglegout)
        newiy = ([l for l in iy if l ∉ danglegout]...,)
        pushfirst!(ycode, EinCode{(newiy,), iy}())
    end
    for i in 1:NT
        ix = ixs[i]
        if !allunique.(ix)
            newix = tunique(ix)
            newxs[i] = einsum(EinCode{(ix,), newix}(), xs[i], size_dict)
            newixs[i] = newix
        end
    end
    if !allunique.(newiy)
        newnewiy = tunique(newiy)
        pushfirst!(ycode, EinCode{(newnewiy,), newiy}())
    end
    res = einsum(EinCode{newixs, newnewiy}(), newxs, size_dict)
    for code in ycode
        # TODO: broadcast and duplicate
        res = einsum(code, (res,), size_dict)
    end
    return res
end

function getalllabels(code::EinCode{ixs, iy}) where {ixs, iy}
    reduce(∪, (ixs..., iy))
end

struct IndexReduction<:OMEinsum.EinRule{1} end

function OMEinsum.match_rule(::IndexReduction, code::EinCode{ixs, iy}) where {ixs, iy}
    length(ixs) > 1 && return false
    ix = ixs[1]
    (OMEinsum.allunique(ix) || !OMEinsum.allunique(iy)) && return false
    allin(iy, ix) && allin(ix, iy) || return false
    return true
end

export allin
function allin(x, y)
    all(xi->xi in y, x)
end

function _get_reductions(ix, iy)
    reds = []
    for l in iy
        count(==(l), ix) > 1 && push!(reds, findall(==(l), ix))
    end
    return reds
end

function OMEinsum.einsum(sm::IndexReduction, code::EinCode{ixs, iy}, xs::Tuple{<:BinarySparseTensor}, size_dict) where {ixs, iy}
    reduce_indices(xs[1], _get_reductions(ixs[1], iy))
end

export allsame
allsame(x::T, mask::T) where T<:Integer = allone(x, mask) || !anyone(x, mask)

function _ymask_from_reds(::Type{Ti}, ndim::Int, reds) where Ti
    ymask = flip(Ti(0), bmask(Ti, 1:ndim))
    for red in reds
        ymask = flip(ymask, bmask(Ti,red[2:end]...))
    end
    return ymask
end

function reduce_indices(t::BinarySparseTensor{Tv,Ti}, masky::Ti, reds) where {Tv,Ti}
    inds = []
    ymask = _ymask_from_reds(Ti, ndims(t), reds)
    bits = baddrs(ymask)
    #masks, ks = BitBasis.group_shift!(ndims(t), position=baddrs(ymask))
    NS = length(ks)
    red_masks = [bmask(Ti, reds...) for red in reds]
    for (ind, val) in zip(t.data.nzind, t.data.nzval)
        b = ind-1
        if all(red->allsame(b, red), red_masks)
            for s in 1:NS
                #@inbounds b = lmove(b, masks[s], ks[s])
                b = readbit(b, bits...)
            end
            push!(inds, b+1)
        end
    end
    return reds
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

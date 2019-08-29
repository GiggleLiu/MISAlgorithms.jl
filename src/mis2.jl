export mis2

"""
Solving MIS problem with sofisticated branching algorithm.
"""
function mis2(eg::EliminateGraph, level::Int=0)
    if nv(eg) == 0
        @show 1
        return 0
    else
        @show 1
        vmin, degmin = mindegree_vertex(eg)
        if degmin <= 1  # DONE
            @show 1.1
            return 1 + eliminate(mis2, eg, NeighborCover(vmin))
        elseif degmin == 2
            @show 1.2
            a, b = Neighbors(eg, vmin)
            if isconnected(eg, a, b)
                @show "1.2.1"
                return 1 + eliminate(mis2, eg, NeighborCover(vmin))
            else
                @show "1.2.2"
                sn = SecondNeighbors(eg, vmin)
                if length(sn) == 1
                    @show "1.2.2.1"
                    w = sn[1]
                    return max(2+eliminate(mis2, eg, NeighborCover(w) ∪ SecondNeighborCover(v)),
                                2+eliminate(mis2, eg, SecondNeighborCover(v)),
                    )
                elseif length(sn) > 1
                    @show "1.2.2.2"
                    return max(eliminate(mis2, eg, NeighborCover(vmin)),
                                eliminate(mis2, eg, MirrorCover(vmin)))
                end
                # NOTE: where is length(sn) == 0?
            end
        elseif degmin == 3 # DONE
            @show 1.3
            a, b, c = Neighbors(eg, vmin)
            nedge = isconnected(eg, a, b) + isconnected(eg, a, c) + isconnected(eg, b, c)
            if nedge == 0
                @show "1.3.1"
                if hasmirror(eg, vmin)
                    @show "1.3.1.1"
                    return max(1+eliminate(mis2, eg, NeighborCover(vmin)),
                                eliminate(mis2, eg, MirrorCover(vmin)))
                else
                    @show "1.3.1.2"
                    return max(1+eliminate(mis2, eg, NeighborCover(vmin)),
                                2 + eliminate(mis2, eg, MirrorCover(a, b)),
                                2 + eliminate(mis2, eg, MirrorCover(a, c) ∪ Vertex(b)),
                                2 + eliminate(mis2, eg, MirrorCover(b, c) ∪ Vertex(a)),
                                )
                end
            elseif nedge == 3
                @show "1.3.2"
                return 1 + eliminate(mis2, eg, NeighborCover(vmin))
            else
                @show "1.3.3"
                return max(1 + eliminate(mis2, eg, NeighborCover(vmin)),
                            eliminate(mis2, eg, MirrorCover(vmin)))
            end
        else # DONE
            @show 1.4
            vmax, degmax = maxdegree_vertex(eg)
            if degmax >= 6 # DONE
                @show "1.4.1"
                return max(1+eliminate(mis2, NeighborCover(vmax)),
                            eliminate(mis2, eg, vmax))
            elseif !isconnected(eg) # DONE
                @show "1.4.2"
                sg = SubGraphWith(eg, v)
                return eliminate(mis2, eg, sg) + mis2(sg)
            elseif dmin == dmax  # DONE
                @show "1.4.3"
                return max(1+eliminate(mis2, eg, NeighborCover(vmax)),
                            eliminate(mis2, eg, MirrorCover(vmax)))
            else
                @show "1.4.4"
                v4, v5 = adjacent45(eg)
                return max(1+eliminate(mis2, eg, NeighborCover(v5)),
                            1+eliminate(mis2, eg, MirrorCover(v5) ∪ NeighborCover(v4)),
                            eliminate(mis2, eg, MirrorCover(v5) ∪ Vertex(v4))
                            )
            end
        end
    end
end

export mis2

"""
Solving MIS problem with sofisticated branching algorithm.
"""
function mis2(graph::EliminateGraph, level::Int=0)
    if nv(graph) == 0
        return 0
    else
        vmin, degmin = mindegree_vertex(eg)
        if degmin <= 1  # DONE
            return 1 + eliminate(mis2, eg, NeighborCover(vmin))
        elseif degmin == 2
            a, b = Neighbors(eg, vmin)
            if isconnected(eg, a, b)
                return 1 + eliminate(mis2, eg, NeighborCover(vmin))
            else
                sn = SecondNeighbors(eg, vmin)
                if length(sn) == 1
                    w = sn[1]
                    return max(2+eliminate(mis2, eg, NeighborCover(w) ∪ SecondNeighborCover(v)),
                                2+eliminate(mis2, eg, SecondNeighborCover(v)),
                    )
                elseif length(sn) > 1
                    return max(eliminate(mis2, eg, NeighborCover(vmin)),
                                eliminate(mis2, eg, MirrorCover(vmin)))
                end
                # NOTE: where is length(sn) == 0?
            end
        elseif degmin == 3 # DONE
            a, b, c = Neighbors(eg, vmin)
            nedge = isconnected(eg, a, b) + isconnected(eg, a, c) + isconnected(eg, b, c)
            if nedge == 0
                if hasmirror(eg, vmin)
                    return max(1+eliminate(mis2, eg, NeighborCover(vmin)),
                                eliminate(mis2, eg, MirrorCover(vmin)))
                else
                    return max(1+eliminate(mis2, eg, NeighborCover(vmin)),
                                2 + eliminate(mis2, eg, MirrorCover(a, b)),
                                2 + eliminate(mis2, eg, MirrorCover(a, c) ∪ Vertex(b)),
                                2 + eliminate(mis2, eg, MirrorCover(b, c) ∪ Vertex(a)),
                                )
                end
            elseif nedge == 3
                return 1 + eliminate(mis2, eg, NeighborCover(vmin))
            else
                return max(1 + eliminate(mis2, eg, NeighborCover(vmin)),
                            eliminate(mis2, eg, MirrorCover(vmin)))
            end
        else # DONE
            vmax, degmax = maxdegree_vertex(eg)
            if degmax >= 6 # DONE
                return max(1+eliminate(mis2, NeighborCover(vmax)),
                            eliminate(mis2, eg, vmax))
            elseif !isconnected(eg) # DONE
                sg = SubGraphWith(eg, v)
                return eliminate(mis2, eg, sg) + mis2(sg)
            elseif dmin == dmax  # DONE
                return max(1+eliminate(mis2, eg, NeighborCover(vmax)),
                            eliminate(mis2, eg, MirrorCover(vmax)))
            else
                v4, v5 = adjacent45(eg)
                return max(1+eliminate(mis2, eg, NeighborCover(v5)),
                            1+eliminate(mis2, eg, MirrorCover(v5) ∪ NeighborCover(v4)),
                            eliminate(mis2, eg, MirrorCover(v5) ∪ Vertex(v4))
                            )
            end
        end
    end
end

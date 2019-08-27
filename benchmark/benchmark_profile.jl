include("mis.jl")

# mask-copy version, 1.5 second
using Random
Random.seed!(2)
eg = rand_eg(60, 0.05)
using BenchmarkTools
#@benchmark neighborcover($eg, 2) seconds=1
@time mis1(eg)
@code_warntype mis1(eg)
@benchmark mis1($eg) seconds=1


using Profile
Profile.clear()
Profile.init(delay=1e-3)
@profile neighborcover(eg, 2)

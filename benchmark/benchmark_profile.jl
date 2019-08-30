using MISAlgorithms

# mask-copy version, 1.5 second
using Random
Random.seed!(2)
eg = rand_eg(60, 0.1)
using BenchmarkTools
#@benchmark neighborcover($eg, 2) seconds=1
@time mis1(eg)
@code_warntype mis1(eg)
@benchmark mis1($eg) seconds=1


using Profile
Profile.clear()
Profile.init(delay=1e-3)
@profile neighborcover(eg, 2)

using Random
Random.seed!(2)
eg = rand_eg(120, 0.1)
@time mis2(eg)

using Profile
eg = rand_eg(80, 0.1)
Profile.clear()
Profile.init(delay=1e-3)
@profile for i=1:30 mis2(eg) end

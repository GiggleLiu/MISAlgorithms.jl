# MISAlgorithms

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://GiggleLiu.github.io/MISAlgorithms.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://GiggleLiu.github.io/MISAlgorithms.jl/dev)
[![Build Status](https://travis-ci.com/GiggleLiu/MISAlgorithms.jl.svg?branch=master)](https://travis-ci.com/GiggleLiu/MISAlgorithms.jl)
[![Codecov](https://codecov.io/gh/GiggleLiu/MISAlgorithms.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/GiggleLiu/MISAlgorithms.jl)

Maximum independent set algorithms, e.g. branching, measure and conqure.

## To develop
Type `]` in a Julia REPL, then input
```julia
pkg> dev git@github.com:GiggleLiu/MISAlgorithms.jl.git
```

## To run an example
```julia
julia> using MISAlgorithms

julia> eg = rand_egraph(60, 0.05);

julia> mis1(eg)  # naive branching algorithm with O(1.44) complexity.
```

## References
* [Exact Exponential Algorithms]()

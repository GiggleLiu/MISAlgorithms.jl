# MISAlgorithms

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

julia> eg = rand_eg(60, 0.05);

julia> mis1(eg)  # naive branching algorithm with O(1.44) complexity.
```

## References
* [Exact Exponential Algorithms](http://www.ii.uib.no/~fomin/BookEA/BookEA.pdf)
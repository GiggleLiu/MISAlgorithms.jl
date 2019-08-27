using Documenter, MISAlgorithms

makedocs(;
    modules=[MISAlgorithms],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/GiggleLiu/MISAlgorithms.jl/blob/{commit}{path}#L{line}",
    sitename="MISAlgorithms.jl",
    authors="JinGuo Liu",
    assets=String[],
)

deploydocs(;
    repo="github.com/GiggleLiu/MISAlgorithms.jl",
)

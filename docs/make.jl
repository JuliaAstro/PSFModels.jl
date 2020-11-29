using PSFKernels
using Documenter

makedocs(;
    modules=[PSFKernels],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaAstro/PSFKernels.jl/blob/{commit}{path}#L{line}",
    sitename="PSFKernels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaAstro.github.io/PSFKernels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/PSFKernels.jl",
    push_preview=true
)

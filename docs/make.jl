using PSFKernels
using Documenter

makedocs(;
    modules=[PSFKernels],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/mileslucas/PSFKernels.jl/blob/{commit}{path}#L{line}",
    sitename="PSFKernels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mileslucas.github.io/PSFKernels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mileslucas/PSFKernels.jl",
)

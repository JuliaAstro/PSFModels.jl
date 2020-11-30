using PSFKernels
using Documenter


setup = quote
    using PSFKernels
end

DocMeta.setdocmeta!(PSFKernels, :DocTestSetup, setup; recursive = true)

makedocs(;
    modules=[PSFKernels],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaAstro/PSFKernels.jl/blob/{commit}{path}#L{line}",
    sitename="PSFKernels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/PSFKernels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API/Reference" => "api.md",
        "Examples" => "examples.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/PSFKernels.jl",
    push_preview=true
)

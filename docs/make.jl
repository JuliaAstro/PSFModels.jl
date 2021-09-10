using PSFModels
using Documenter


setup = quote
    using PSFModels
end

DocMeta.setdocmeta!(PSFModels, :DocTestSetup, setup; recursive = true)

makedocs(;
    modules=[PSFModels],
    authors="Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo="https://github.com/JuliaAstro/PSFModels.jl/blob/{commit}{path}#L{line}",
    sitename="PSFModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaastro.github.io/PSFModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API/Reference" => "api.md",
        "Examples" => "examples.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaAstro/PSFModels.jl",
    push_preview=true,
    devbranch="main"
)

using PSFModels
using Documenter
using Documenter.Remotes: GitHub

setup = quote
    using PSFModels
end

DocMeta.setdocmeta!(PSFModels, :DocTestSetup, setup; recursive = true)
include("pages.jl")
makedocs(;
    modules = [PSFModels],
    authors = "Miles Lucas <mdlucas@hawaii.edu> and contributors",
    repo = GitHub("JuliaAstro/PSFModels.jl"),
    sitename = "PSFModels.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://juliaastro.github.io/PSFModels.jl",
        assets = String[],
        canonical="https://juliaastro.org/PSFModels/stable/",
    ),
    pages = pages,
)

deploydocs(;
    repo = "github.com/JuliaAstro/PSFModels.jl",
    push_preview = true,
    devbranch = "main",
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)

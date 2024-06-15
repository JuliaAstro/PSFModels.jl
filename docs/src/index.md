```@meta
CurrentModule = PSFModels
```

# PSFModels.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/PSFModels.jl)
[![Build Status](https://github.com/juliaastro/PSFModels.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliaastro/PSFModels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFModels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliaastro/PSFModels.jl)
[![License](https://img.shields.io/github/license/juliaastro/PSFModels.jl?color=yellow)](https://github.com/juliaastro/PSFModels.jl/blob/main/LICENSE)

## Installation

PSFModels can be added from the Julia package manager

```julia-repl
julia>]

(@v1.6) pkg> add PSFModels
```

## Getting Started

To import the library

```julia-repl
julia> using PSFModels
```

None of the models are exported to avoid namespace clashes, but it can be verbose to continuously rewrite `PSFModels`. You can either import names directly

```julia-repl
julia> using PSFModels: gaussian

julia> model = gaussian(x=0, y=0, fwhm=8)
```

or you can create an alias for `PSFModels`

```julia
# julia version 1.5 or below
using PSFModels
const M = PSFModels
# julia version 1.6 or above
import PSFModels as M

model = M.gaussian(x=0, y=0, fwhm=10)
```

```@docs
PSFModels
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/PSFModels.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/PSFModels.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaAstro/PSFModels.jl/issues).

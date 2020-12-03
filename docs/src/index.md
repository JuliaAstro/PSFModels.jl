```@meta
CurrentModule = PSFModels
```

# PSFModels.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/PSFModels.jl)
[![Build Status](https://github.com/juliaastro/PSFModels.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/PSFModels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFModels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFModels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/PSFModels.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

PSFModels can be added from the Julia package manager

```julia
julia>]

(@v1.6) pkg> add PSFModels
```

## Getting Started

To import the library

```julia
julia> using PSFModels
```

None of the kernels are exported to avoid namespace clashes, but it can be verbose. You can either import names directly

```julia
julia> using PSFModels: Gaussian

julia> kernel = Gaussian(8)
```

or you can create an alias for `PSFModels`

```julia
# julia version 1.5 or below
using PSFModels
const kerns = PSFModels
# julia version 1.6 or above
using PSFModels as kerns

kernel = kerns.Gaussian(10)
```

```@docs
PSFModels
```

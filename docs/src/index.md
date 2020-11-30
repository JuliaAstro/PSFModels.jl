```@meta
CurrentModule = PSFKernels
```

# PSFKernels.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/PSFKernels.jl)
[![Build Status](https://github.com/juliaastro/PSFKernels.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/PSFKernels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFKernels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFKernels.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/PSFKernels.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

PSFKernels can be added from the Julia package manager

```julia
julia>]

(@v1.6) pkg> add PSFKernels
```

to import the library

```julia
julia> using PSFKernels
```

None of the kernels are exported to avoid namespace clashes, but it can be verbose. You can either import names directly

```julia
julia> using PSFKernels: Gaussian

julia> kernel = Gaussian(8)
```

or you can create an alias for `PSFKernels`

```julia
# julia version 1.5 or below
using PSFKernels
const k = PSFKernels
# julia version 1.6 or above
using PSFKernels as k

kernel = k.Gaussian(10)
```

```@docs
PSFKernels
```

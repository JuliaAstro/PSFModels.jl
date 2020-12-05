# PSFModels.jl

[![Build Status](https://github.com/juliaastro/PSFModels.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/PSFModels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFModels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFModels.jl/branch/master/graph/badge.svg?branch=master)](https://codecov.io/gh/juliaastro/PSFModels.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/PSFModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/PSFModels.jl/dev)

Fast, allocation-free point-spread function (PSF) representations

## Models

* `PSFModels.Gaussian` (or `Normal`)
* `PSFModels.AiryDisk`
* `PSFModels.Moffat`

## Installation

From the Julia REPL

```julia
julia> ]

(@v1.5) pkg> add PSFModels
```

To import the library

```julia
julia> using PSFModels
```

None of the models are exported to avoid namespace clashes, but it can be verbose. You can either import names directly

```julia
julia> using PSFModels: Gaussian

julia> model = Gaussian(8)
```

or you can create an alias for `PSFModels`

```julia
# julia version 1.5 or below
using PSFModels
const M = PSFModels
# julia version 1.6 or above
using PSFModels as M

model = M.Gaussian(10)
```

## Usage

For more in-depth usage and examples, please see the [documentation](https://juliaastro.github.io/PSFModels.jl/dev/).

```julia
using PSFModels

m = PSFModels.Gaussian(8)           # bivariate gaussian with a FWHM of 8 pixels
m = PSFModels.Gaussian((7.4, 8.2))  # specify FWHM for each axis
m = PSFModels.Gaussian([1 0; 0 1])  # specify FWHM as a correlated matrix

m = PSFModels.Gaussian(12, 25, 8.2) # specifiy location in pixel coordinates
m = PSFModels.Gaussian([12, 25], 8.2)

mf0 = PSFModels.Gaussian{Float32}(8.2) # output guaranteed to be Float32
```

```julia
m[0, 0]      # "index" the model at [y, x]
m[:, 0]
m(0.3, 1.0)  # directly query value at (x, y)
m([1.2, 0.4])

# scalar multiplication or division will create a ScaledPSFModel
20 * m # or `m * 20`
m / 20

# evaluate `m` over its indices forming an array
collect(m)

# broadcasting will take the axes of the other arrays
arr = randn(101, 101)
m .* arr

## (nearly) allocation-free loss function
# get overlapped cutouts for the PSF and the array
inds = map(intersect, axes(arr), axes(m))
arr_stamp = @view arr[inds...]
# use `amp * m` to create ScaledPSFModel (could also just broadcast)
amp = 1.24
m_stamp = @view (amp * m)[inds...]
resid = sum(abs2, arr_stamp .- m_stamp) # chi-square loss
```

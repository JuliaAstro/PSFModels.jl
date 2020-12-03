# PSFModels.jl

[![Build Status](https://github.com/juliaastro/PSFModels.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/PSFModels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFModels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFModels.jl/branch/master/graph/badge.svg?branch=master)](https://codecov.io/gh/juliaastro/PSFModels.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/PSFModels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/PSFModels.jl/dev)

Fast, allocation-free point-spread function (PSF) representations

## Kernels

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

## Usage

For more in-depth usage and examples, please see the [documentation](https://juliaastro.github.io/PSFModels.jl/dev/).

```julia
using PSFModels

k = PSFModels.Gaussian(8)           # bivariate gaussian with a FWHM of 8 pixels
k = PSFModels.Gaussian((7.4, 8.2))  # specify FWHM for each axis
k = PSFModels.Gaussian([1 0; 0 1])  # specify FWHM as a correlated matrix

k = PSFModels.Gaussian(12, 25, 8.2) # specifiy location in pixel coordinates
k = PSFModels.Gaussian([12, 25], 8.2)

kf0 = PSFModels.Gaussian{Float32}(8.2) # output guaranteed to be Float32
```

```julia
k[0, 0]      # "index" the kernel
k[:, 0]
k(0.3, 1.0)  # directly query value
k([1.2, 0.4])

# evaluate `k` over its indices forming an array
collect(k)

# broadcasting will take the axes of the other arrays
arr = randn(101, 101)
k .* arr

## (nearly) allocation-free loss function
# get overlapped cutouts for the PSF and the array
inds = map(intersect, axes(arr), axes(k))
arr_stamp = @view arr[inds...]
kern_stamp = @view k[inds...]
amp = 1.24
resid = sum(abs2, arr_stamp .- amp .* kern_stamp) # chi-square loss
```

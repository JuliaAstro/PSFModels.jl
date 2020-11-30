# PSFKernels.jl

[![Build Status](https://github.com/juliaastro/PSFKernels.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/PSFKernels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFKernels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFKernels.jl/branch/master/graph/badge.svg?branch=master)](https://codecov.io/gh/juliaastro/PSFKernels.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/PSFKernels.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/PSFKernels.jl/dev)

Fast, allocation-free point-spread function (PSF) representations

## Kernels

* `PSFKernels.Gaussian` (or `Normal`)
* `PSFKernels.AiryDisk`
* `PSFKernels.Moffat`

## Installation

From the Julia REPL

```julia
julia> ]

(@v1.5) pkg> add PSFKernels
```

## Usage

For more in-depth usage and examples, please see the [documentation]([https://](https://juliaastro.github.io/PSFKernels.jl/dev)).

```julia
using PSFKernels

k = PSFKernels.Gaussian(8)           # bivariate gaussian with a FWHM of 8 pixels
k = PSFKernels.Gaussian((7.4, 8.2))  # specify FWHM for each axis
k = PSFKernels.Gaussian([1 0; 0 1])  # specify FWHM as a correlated matrix

k = PSFKernels.Gaussian(12, 25, 8.2) # specifiy location in pixel coordinates
k = PSFKernels.Gaussian([12, 25], 8.2)

kf0 = PSFKernels.Gaussian{Float32}(8.2) # output guaranteed to be Float32
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

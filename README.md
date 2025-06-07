# PSFModels.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.org/PSFModels/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.org/PSFModels.jl/dev/)

[![CI](https://github.com/JuliaAstro/PSFModels.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaAstro/PSFModels.jl/actions/workflows/ci.yml)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFModels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/juliaastro/PSFModels.jl/graph/badge.svg?token=Jy06ZuwvVi)](https://codecov.io/gh/juliaastro/PSFModels.jl)
![License](https://img.shields.io/github/license/JuliaAstro/PSFModels.jl?color=yellow)

Fast, allocation-free point-spread function (PSF) representations

## Models

* `gaussian` (or `normal`)
* `airydisk`
* `moffat`

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

```julia
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

model = M.gaussian(fwhm=10)
```

## Usage

For more in-depth usage and examples, please see the [documentation](https://juliaastro.github.io/PSFModels.jl/dev/).

First, load the package

```julia
julia> using PSFModels
```

### Evaluating models

Directly evaluating the functions is the most straightforward way to use this package

```julia
julia> gaussian(0, 0; x=0, y=0, fwhm=3)
1.0

julia> gaussian(BigFloat, 0, 0; x=0, y=0, fwhm=3, amp=0.1)
0.1000000000000000055511151231257827021181583404541015625
```

We also provide "curried" versions of the functions, which allow you to specify the parameters and evaluate the PSF later

```julia
julia> model = gaussian(x=0, y=0, fwhm=3);

julia> model(0, 0)
1.0
```

If we want to collect the model into a dense matrix, simply iterate over indices

```julia
julia> inds = CartesianIndices((-2:2, -2:2));

julia> model.(inds) # broadcasting
5×5 Matrix{Float64}:
 0.0850494  0.214311  0.291632  0.214311  0.0850494
 0.214311   0.54003   0.734867  0.54003   0.214311
 0.291632   0.734867  1.0       0.734867  0.291632
 0.214311   0.54003   0.734867  0.54003   0.214311
 0.0850494  0.214311  0.291632  0.214311  0.0850494
```

This makes it very easy to evaluate the PSF on the same axes as an image (array)

```julia
julia> img = randn(5, 5);

julia> model.(CartesianIndices(img))
5×5 Matrix{Float64}:
 0.54003      0.214311     0.0459292    0.00531559   0.000332224
 0.214311     0.0850494    0.018227     0.00210949   0.000131843
 0.0459292    0.018227     0.00390625   0.000452087  2.82555e-5
 0.00531559   0.00210949   0.000452087  5.2322e-5    3.27013e-6
 0.000332224  0.000131843  2.82555e-5   3.27013e-6   2.04383e-7
```

this is trivially expanded to fit "stamps" in images

```julia
julia> big_img = randn(1000, 1000);

julia> stamp_inds = (750:830, 400:485);

julia> stamp = @view big_img[stamp_inds...];

julia> stamp_model = model.(CartesianIndices(stamp_inds));
```

or we can create a loss function for fitting PSFs without allocating any memory. We are simply iterating over the image array!

```julia
julia> using Statistics

julia> mse = mean(I -> (big_img[I] - model(I))^2, CartesianIndices(stamp_inds));
```

### Fitting data

There exists a simple, yet powerful, API for fitting data with these PSF models. See the [full documentation](https://juliaastro.github.io/PSFModels.jl/dev) for more details and examples.

```julia
# `fit` is not exported to avoid namespace clashes
using PSFModels: fit

data = # load data
stamp_inds = # optionally choose indices to "cutout"

# use an isotropic Gaussian
P0 = (x=12, y=13, fwhm=3.2, amp=0.1)
params, synthpsf = fit(gaussian, P0, data, stamp_inds)

# elliptical, rotated Gaussian
P0 = (x=12, y=13, fwhm=(3.2, 3.2), amp=0.1, theta=0)
params, synthpsf = fit(gaussian, P0, data, stamp_inds)

# obscured Airy disk
P0 = (x=12, y=13, fwhm=3.2, amp=0.1, ratio=0.3)
params, synthpsf = fit(airydisk, P0, data, stamp_inds)

# bivariate Moffat with arbitrary alpha
P0 = (x=12, y=13, fwhm=(3.2, 3.2), amp=0.1, alpha=1)
# fixed ("frozen") rotation angle
func_kwargs = (;theta=15)
params, synthpsf = fit(moffat, P0, data, stamp_inds; func_kwargs)
```

### Plotting models

We provide simple user recipes from [RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl), which can be called with `psfplot`/`psfplot!`

```julia
using Plots

inds = (1:30, 1:30)
model = airydisk(x=12, y=13, fwhm=(4.5, 6.7), theta=12, ratio=0.3)
psfplot(model, inds, colorbar_scale=:log10)
```

## Contributing and Support

If you would like to contribute, feel free to open a [pull request](https://github.com/JuliaAstro/PSFModels.jl/pulls). If you want to discuss something before contributing, head over to [discussions](https://github.com/JuliaAstro/PSFModels.jl/discussions) and join or open a new topic. If you're having problems with something, please open an [issue](https://github.com/JuliaAstro/PSFModels.jl/issues).

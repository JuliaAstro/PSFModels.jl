```@meta
DocTestSetup = quote
    using PSFModels
end
```

# Introduction

Statistical models for constructing point-spread functions (PSFs).

## Models

The following models are currently implemented
* [`gaussian`](@ref)/[`normal`](@ref)
* [`airydisk`](@ref)
* [`moffat`](@ref)

## Parameters

In general, the PSFs have a position, a full-width at half-maximum (FWHM)
measure, and an amplitude. The position follows a 1-based pixel coordinate
system, where `(1, 1)` represents the *center* of the bottom left pixel. This
matches the indexing style of Julia as well as DS9, IRAF, SourceExtractor, and
WCS. The FWHM is a consistent scale parameter for the models. That means a
`gaussian` with a FWHM of 5 will be visually similar to an `airydisk` with a
FWHM of 5. All models support a scalar (isotropic) FWHM and a FWHM for each axis
(diagonal), as well as arbitrarily rotating the PSF.

!!! warning "Pixel-convention"

    The pixel convention adopted here is that the bottom-left pixel *center* is
    `(1, 1)`. The column-major memory layout of julia puts the `x` axis as the
    rows of a matrix and the `y` axis as the columns. In other words, the axes
    unpack like
    ```julia
    xs, ys = axes(image)
    ```


## Usage

### Evaluating models

Directly evaluating the functions is the most straightforward way to use this package

```jldoctest model
julia> gaussian(0, 0; x=0, y=0, fwhm=3)
1.0

julia> gaussian(BigFloat, 0, 0; x=0, y=0, fwhm=3, amp=0.1, bkg=1)
1.100000000000000088817841970012523233890533447265625
```

We also provide "curried" versions of the functions, which allow you to specify
the parameters and evaluate the PSF later

```jldoctest model
julia> model = gaussian(x=0, y=0, fwhm=3);

julia> model(0, 0)
1.0
```

If we want to collect the model into a dense matrix, simply iterate over indices

```jldoctest model
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

```jldoctest model
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

```jldoctest model
julia> big_img = randn(1000, 1000);

julia> stamp_inds = (750:830, 400:485);

julia> stamp = @view big_img[stamp_inds...];

julia> stamp_model = model.(CartesianIndices(stamp_inds));
```

or we can create a loss function for fitting PSFs without allocating any memory.
We are simply iterating over the image array!

```jldoctest model
julia> using Statistics

julia> mse = mean(I -> (big_img[I] - model(I))^2, CartesianIndices(stamp_inds));
```

### Fitting data

There exists a simple, yet powerful, API for fitting data with
[`PSFModels.fit`](@ref).

```julia
# `fit` is not exported to avoid namespace clashes
using PSFModels: fit

data = # load data
stamp_inds = # optionally choose indices to "cutout"

# use an isotropic Gaussian
params, synthpsf = fit(gaussian, (x=12, y=13, fwhm=3.2, amp=0.1),
                       data, stamp_inds)

# elliptical, rotated Gaussian
params, synthpsf = fit(gaussian, (x=12, y=13, fwhm=(3.2, 3.2), amp=0.1, theta=0)
                       data, stamp_inds)

# obscured Airy disk
params, synthpsf = fit(airydisk, (x=12, y=13, fwhm=3.2, amp=0.1, ratio=0.3),
                       data, stamp_inds)

# bivariate Moffat with arbitrary alpha
params, synthpsf = fit(moffat, (x=12, y=13, fwhm=(3.2, 3.2), amp=0.1, alpha=1),
                       data, stamp_inds)
```

### Plotting

Finally, we provide plotting recipes (`psfplot`/`psfplot!`) from
[RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl), which can be
seen in use in the [API/Reference](@ref) section.

```julia
using Plots

model = gaussian(x=0, y=0, fwhm=(8, 10), theta=12)
psfplot(model, -30:30, -30:30, colorbar_scale=:log10)
```

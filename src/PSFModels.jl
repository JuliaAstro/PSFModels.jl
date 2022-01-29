"""
# PSFModels

Statistical models for constructing point-spread functions (PSFs). These models act like matrices but without allocating any memory, which makes them efficient to fit and apply.

## Models

The following models are currently implemented
* [`PSFModels.Gaussian`](@ref)/[`PSFModels.Normal`](@ref)
* [`PSFModels.AiryDisk`](@ref)
* [`PSFModels.Moffat`](@ref)

## Parameters

In general, the PSFs have a position, a full-width at half-maximum (FWHM) measure, and an amplitude. The position follows a 1-based pixel coordinate system, where `(1, 1)` represents the *center* of the bottom left pixel. This matches the indexing style of Julia as well as DS9, IRAF, SourceExtractor, and WCS. If a position is not specified, it is set to `(0, 0)`. The FWHM is a consistent scale parameter for the models. All models support a scalar (isotropic) FWHM and a FWHM for each axis (diagonal).

## Usage

Using the models should feel just like an array. In fact, `PSFModels.PSFModel <: AbstractMatrix`. However, no data is stored and no allocations have to be made. In other words, representing the models as matrices is merely a convenience, since typically astronomical data is stored in dense arrays. Another way of thinking of these is a lazy array that applies a function when indexed, rather than returning data stored in memory.

```jldoctest model
julia> m = PSFModels.Gaussian(fwhm=3); # centered at (0, 0)


julia> m[0, 0] # [x, y] for indexing
1.0

julia> m(0, 0) # (x, y) for evaluating
1.0
```

Because the model is a matrix, it needs to have a size. Each model has a bounding box which can be controlled with the `extent` keyword. By default the extent is set by a scalar factor of the FWHM (e.g., `maxsize * FWHM` pixels), centered around the PSF, and rounded up. We can see how this alters the indices from a typical `Matrix`

```jldoctest model
julia> size(m) # default 'stamp' size is fwhm * 3
(11, 11)

julia> axes(m)
(-5:5, -5:5)
```

if we want to collect the model into a dense matrix, regardless of the indexing (e.g. to prepare for cross-correlation), we can simply

```jldoctest model
julia> stamp = collect(m);

```

these axes are merely a convenience for bounding the model, since they accept any real number as input.

```jldoctest model
julia> m[100, 10000] # index-like inputs [x, y]
0.0

julia> m(2.4, 1.7) # valid for any real (x, y)
0.0696156536973086
```

By bounding the model, we get a cutout which can be applied to arrays with much larger dimensions without having to iterate over the whole matrix

```jldoctest
julia> big_mat = ones(1001, 1001);

julia> model = PSFModels.Gaussian(x=51, y=51, fwhm=2);


julia> ax = map(intersect, axes(big_mat), axes(model))
(48:54, 48:54)

julia> cutout = @view big_mat[ax...]
7×7 view(::Matrix{Float64}, 48:54, 48:54) with eltype Float64:
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> stamp = @view model[ax...];


julia> photsum = sum(cutout .* stamp)
4.5322418212890625
```
Nice- we only had to reduce ~50 pixels instead of ~1,000,000 to calculate the aperture sum, all in under a microsecond (on my machine).

Since the models are lazy, that means the type of the output can be specified, as long as it can be converted to from a real number (so no integer types).

```jldoctest
julia> mbig = PSFModels.Gaussian(BigFloat, fwhm=12);


julia> sum(mbig)
163.07467408408593885562554918859656805096847165259532630443572998046875
```

finally, we provide plotting recipes from [RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl), which can be seen in use in the [API/Reference](@ref) section.

```julia
using Plots
model = PSFModels.Gaussian(8)
plot(model)              # default axes
plot(model, 1:5, :)    # custom axes (x, y)
plot(model, axes(other)) # use axes from other array
```
"""
module PSFModels

using CoordinateTransformations
using Distances
using KeywordCalls
using Rotations
using SpecialFunctions
using StaticArrays

include("core.jl")
include("functions/gaussian.jl")
include("functions/airy.jl")
# include("functions/moffat.jl")
# include("plotting.jl")

end # module PSFModels

"""
# PSFModels

Statistical models for constructing point-spread functions (PSFs). These models act like matrices but without allocating any memory, which makes them efficient to fit and apply.

## Models

The following models are currently implemented
* [`PSFModels.Gaussian`](@ref)
* [`PSFModels.AiryDisk`](@ref)
* [`PSFModels.Moffat`](@ref)

## Parameters

In general, the PSFs have a position, a full-width at half-maximum (FWHM) measure, and an amplitude. The position a 1-based pixel coordinate system, where `(1, 1)` represents the *center* of the bottom left pixel. This matches the indexing style of Julia as well as DS9, IRAF, SourceExtractor, and WCS. If a position is not specified, it is set to `(0, 0)`. The FWHM is a consistent scale parameter for the models. All models support a scalar (isotropic) FWHM and a FWHM for each axis (diagonal).

## Usage

Using the models should feel just like an array. In fact, `PSFModels.PSFModel <: AbstractMatrix`. However, no data is stored and no allocations have to be made. In other words, representing the models as matrices is merely a convenience, since typically astronomical data is stored in dense arrays.

```jldoctest model
julia> m = PSFModels.Gaussian(5); # fwhm of 5 pixels, centered at (0, 0)

julia> m[0, 0] # [y, x] for indexing
1.0

julia> m(0, 0) # (x, y) for evaluating
1.0
```

!!! note "axis order"

    It's important to note the difference in the axis ordering between the index-style calls and the function-style calls. The index-style calls are reverse cartesian order (e.g., `(z, y, x)`), while function calls are the typical cartesian order `(x, y, z)`. Regardless, the constructors are always in cartesian order (`(x, y, z)`).

To control the amplitude, the best method is using scalar multiplication or division. These operations create another lazy object ([`ScaledPSFModel`](@ref)) that scales the original model without having to broadcast and potentially allocate.

```jldoctest model
julia> m_scaled = 20 * m;

julia> m_scaled(0, 0)
20.0

julia> m′ = m_scaled / 20;

julia> m′(0, 0)
1.0
```

Because the model is a matrix, it needs to have a size. In this case, the size is `maxsize * FWHM` pixels, centered around the origin, and rounded up. We can see how this alters the indices from a typical `Matrix`

```jldoctest model
julia> size(m)
(17, 17)

julia> axes(m)
(-8:8, -8:8)
```

if we want to collect the model into a dense matrix, regardless of the indexing (e.g. to prepare for cross-correlation), we can simply

```jldoctest model
julia> stamp = collect(m);
```

these axes are merely a convenience for bounding the model, since they accept any real number as input.

```jldoctest model
julia> m[100, 10000] # index-like inputs [y, x]
0.0

julia> m(2.4, 1.7) # valid for any real (x, y)
0.38315499005194587
```

By bounding the model, we get a cutout which can be applied to arrays with much larger dimensions without having to iterate over the whole matrix

```jldoctest
julia> big_mat = ones(101, 101);

julia> model = PSFModels.Gaussian(51, 51, 2); # center of big_mat, fwhm=2

julia> ax = map(intersect, axes(big_mat), axes(model))
(48:54, 48:54)

julia> cutout = @view big_mat[ax...]
7×7 view(::Array{Float64,2}, 48:54, 48:54) with eltype Float64:
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
Nice- we only had to reduce ~50 pixels instead of ~10,000 to calculate the aperture sum, all in under a microsecond (on my machine).

Since the models are lazy, that means the type of the output can be specified, as long as it can be converted to from a real number (so no integer types).

```jldoctest
julia> mbig = PSFModels.Gaussian{BigFloat}(12);

julia> sum(mbig)
163.07467408408593790971336380361822460116627553361468017101287841796875
```

finally, we provide plotting recipes from [RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl), which can be seen in use in the [API/Reference](@ref) section.

```julia
using Plots
model = PSFModels.Gaussian(8)
plot(model)              # default axes
plot(model, 1:5, 1:5)    # custom axes
plot(model, axes(other)) # use axes from other array
```

!!! tip "Tip: Automatic Differentation"
    Forward-mode AD libraries tend to use dual numbers, which can cause headaches getting the types correct. We recommend using the *primal vector*'s element type to avoid these headaches
    ```julia
    # example generative model for position and scalar fwhm
    model(X::AbstractVector{T}) where {T} = PSFModels.Gaussian{T}(X...)
    ```
"""
module PSFModels

using CoordinateTransformations
using Distances
using SpecialFunctions
using StaticArrays

"""
    PSFModels.PSFModel{T} <: AbstractMatrix{T}

Abstract type for PSF models.

In general, all `PSFModel`s have a set of pre-determined axes (the size is set upon creation) but they are lazy. That is, no memory is allocated and the values are calculated on the fly.

# Interface
The interface to define a model is as follows (for an example model `Model`)

| method | description |
|:-------|:------------|
| `Model()` | constructor(s) |
| `Base.size(m::Model)` | size, necessary for `AbstractArray` interface |
| `Base.axes(m::Model)` | axes, necessary for `AbstractArray` interface |
| `(m::Model)(point::AbstractVector)` | evaluate the model at the point in 2d space (x, y) |

browsing through the implementation of [`PSFModels.Gaussian`](@ref) should give a good idea of how to create a model
"""
abstract type PSFModel{T} <: AbstractMatrix{T} end

# always inbounds
Base.checkbounds(::Type{Bool}, ::PSFModel, idx...) = true
Base.checkbounds(::Type{Bool}, ::PSFModel, idx::CartesianIndex) = true

# in general, parse to static vector
(model::PSFModel)(point...) = model(SVector(point))
(model::PSFModel)(point::Tuple) = model(SVector(point))
# disallow index to avoid confusion 
(model::PSFModel)(::CartesianIndex) = error("PSF models should be indexed using `getindex`, (equivalently `[]`)")

# getindex just calls model with reversed indices
Base.getindex(model::PSFModel, idx::Vararg{<:Number,2}) = model(reverse(idx))

# broadcasting hack to slurp other axes (doesn't work for numbers)
Broadcast.combine_axes(kern::PSFModel, other) = axes(other)
Broadcast.combine_axes(other, kern::PSFModel) = axes(other)

function indices_from_extent(pos, fwhm, maxsize)
    halfextent = @. maxsize * fwhm / 2
    lower = @. floor(Int, pos - halfextent)
    upper = @. ceil(Int, pos + halfextent)
    return (first(lower):first(upper), last(lower):last(upper))
end

function indices_from_extent(pos, extent)
    halfextent = @. extent / 2
    lower = @. round(Int, pos - halfextent)
    upper = @. round(Int, pos + halfextent)
    return (first(lower):first(upper), last(lower):last(upper))
end

# TODO
# function indices_from_extent(pos, fwhm::AbstractMatrix, maxsize)
#     halfextent = maxsize .* fwhm ./ 2
#     lower = @. floor(Int, pos - halfextent)
#     upper = @. ceil(Int, pos - halfextent)
# end

include("gaussian.jl")
include("moffat.jl")
include("airy.jl")
include("scaled.jl")
include("plotting.jl")

end # module PSFModels

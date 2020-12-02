"""
# PSFKernels

Statistical kernels for constructing point-spread functions (PSFs). These kernels act like matrices but without allocating any memory, which makes them efficient to fit and apply.

## Kernels

The following kernels are currently implemented
* [`PSFKernels.Gaussian`](@ref)
* [`PSFKernels.AiryDisk`](@ref)
* [`PSFKernels.Moffat`](@ref)

## Usage

Using the kernels should feel just like an array. In fact, `PSFKernels.PSFKernel <: AbstractMatrix`. However, no data is stored and no allocations have to be made. In other words, representing the kernels as matrices is merely a convenience, since typically astronomical data is stored in dense arrays.

```jldoctest kernel
julia> k = PSFKernels.Gaussian(5); # fwhm of 5 pixels, centered at (0, 0)

julia> k[0, 0]
1.0
```
because the kernel is a matrix, it needs to have a size. In this case, the size is `maxsize * FWHM` pixels, centered around the origin, and rounded up. We can see how this alters the indices from a typical `Matrix`

```jldoctest kernel
julia> size(k)
(17, 17)

julia> axes(k)
(-8:8, -8:8)
```

if we want to collect the kernel into a dense matrix, regardless of the indexing (e.g. to prepare for cross-correlation), we can simply

```jldoctest kernel
julia> stamp = collect(k);
```

these axes are merely a convenience for bounding the kernel, since they accept any real number as input. 

```jldoctest kernel
julia> k[100, 10000] # valid for index-like inputs
0.0

julia> k(2.4, 1.7) # valid for any number
0.38315499005194587
```

By bounding the kernel, we get a cutout which can be applied to arrays with much larger dimensions without having to iterate over the whole matrix

!!! note "Coordinate System"
    The PSFs follows a 1-based image coordinate system, where `(1, 1)` represents the *center* of the bottom left pixel. This corresponds with Julia indexing, as well as DS9 and IRAF.

```jldoctest
julia> big_mat = ones(101, 101);

julia> small_kern = PSFKernels.Gaussian(51, 51, 2); # center of big_mat, fwhm=2

julia> ax = map(intersect, axes(big_mat), axes(small_kern))
(48:54, 48:54)

julia> cutout = big_mat[ax...]
7×7 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> photsum = sum(cutout .* small_kern[ax...])
4.5322418212890625
```
Nice- we only had to reduce ~50 pixels instead of ~10,000 to calculate the aperture sum, and with some care we could make it allocation-free.

Since the kernels are lazy, that means the type of the output can be specified, as long as it can be converted to from a real number (so no integer types).

```jldoctest
julia> kbig = PSFKernels.Gaussian{BigFloat}(12);

julia> sum(kbig)
163.07467408408593790971336380361822460116627553361468017101287841796875
```

finally, we provide plotting recipes from [RecipesBase.jl](https://github.com/JuliaPlots/RecipesBase.jl), which can be seen in use in the [API/Reference](@ref) section.

```julia
using Plots
kernel = PSFKernels.Gaussian(8)
plot(kernel)              # default axes
plot(kernel, 1:5, 1:5)    # custom axes
plot(kernel, axes(other)) # use axes from other array
```

!!! tip "Tip: Automatic Differentation"
    Forward-mode AD libraries tend to use dual numbers, which can cause headaches getting the types correct. We recommend using the primal vector's element type to avoid these headaches
    ```julia
    # example generative model for position and scalar fwhm
    kernel(X::AbstractVector{T}) where {T} = PSFKernels.Gaussian{T}(X...)
    ```
"""
module PSFKernels

using CoordinateTransformations
using Distances
using SpecialFunctions
using StaticArrays

"""
    Kernels.PSFKernel{T} <: AbstractMatrix{T}

Abstract type for PSF Kernels.

In general, all `PSFKernel`s have a set of pre-determined axes (the size is set upon creation) but they are lazy. That is, no memory is allocated and the values are calculated on the fly.

# Interface
The interface to define a kernel is as follows (for an example kernel `Kernel`)

| method | description |
|:-------|:------------|
| `Kernel()` | constructor(s) |
| `Base.size(k::Kernel)` | size, necessary for `AbstractArray` interface |
| `Base.axes(k::Kernel)` | axes, necessary for `AbstractArray` interface |
| `(k::Kernel)(point::AbstractVector)` | evaluate the kernel at the point in 2d space |

browsing through the implementation of [`PSFKernels.Gaussian`](@ref) should give a good idea of how to create a kernel
"""
abstract type PSFKernel{T} <: AbstractMatrix{T} end

# always inbounds
Base.checkbounds(::Type{Bool}, ::PSFKernel, idx...) = true
Base.checkbounds(::Type{Bool}, ::PSFKernel, idx::CartesianIndex) = true

# in general, parse to static vector
(kernel::PSFKernel)(point...) = kernel(SVector(point))
(kernel::PSFKernel)(point::Tuple) = kernel(SVector(point))
(kernel::PSFKernel)(idx::CartesianIndex) = kernel(SVector(idx.I))

# getindex just calls kernel
Base.getindex(kernel::PSFKernel, idx::Vararg{Int,2}) = kernel(idx)

# broadcasting hack to slurp other axes (doesn't work for numbers)
Broadcast.combine_axes(kern::PSFKernel, other) = axes(other)
Broadcast.combine_axes(other, kern::PSFKernel) = axes(other)

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
include("plotting.jl")

end # module PSFKernels

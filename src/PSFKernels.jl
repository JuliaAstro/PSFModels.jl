"""
    PSFKernels

Statistical kernels for constructing point-spread functions (PSFs). These kernels act like matrices but without allocating any memory, which makes them efficient to fit and apply.

# Kernels

The following kernels are currently implemented
* [`PSFKernels.Gaussian`](@ref)
* [`PSFKernels.AiryDisk`](@ref)
* [`PSFKernels.Moffat`](@ref)

# Usage

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
julia> k[100, 10000]
0.0
```

By bounding the kernel, we get a cutout which can be applied to arrays with much larger dimensions without having to iterate over the whole matrix

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

!!! tip "Tip: Automatic Differentation"
Forward-mode AD libraries tend to use dual numbers, which can cause headaches getting the types correct. We recommend using the primal vector's element type to avoid these headaches
```julia
# example generative model for position and scalar fwhm
kernel(X::AbstractVector{T}) where {T} = Kernels.Gaussian{T}(X...)
```

# Examples

Here is a brief example which shows how to construct a loss function for fitting a `PSFKernel` to some data.

```julia
data = ones(101, 101)
# generative model
function kernel(X::AbstractVector{T}) where T
    position = X[1:2] # x, y position
    fwhm = X[3:4]     # fwhm_x, fwhm_y
    return PSFKernels.Gaussian{T}(position, fwhm)
end
# objective function
function loss(X::AbstractVector, data)
    k = kernel(X)
    amp = X[5]
    # l2-distance loss (χ² loss)
    return sum(abs2, data .- amp .* k[axes(data)...])
end

test_params = Float64[51, 51, 10, 10, 1]
loss(test_params, data)

# output

10031.03649468148
```

The objective function can then be used with an optimization library like [Optim.jl](https://github.com/JuliaOpt/Optim.jl) to find best-fitting parameters
```julia
using Optim

# Fit our data using test_params as a starting point
# uses Nelder-Mead optimization
res = optimize(P -> loss(P, data), test_params)

# utilize automatic differentiation (AD) to enable
# advanced algorithms, like Newton's method
res_ad = optimize(P -> loss(P, data), test_params, Newton(); autodiff=:forward)
```
we can see which result has the better loss, and then use the generative model to create a kernel that we can use elsewhere
```julia
best_res = minimum(res) < minimum(res_ad) ? res : res_ad
best_fit_params = Optim.minimizer(best_res)
model = kernel(best_fit_params)
```
"""
module PSFKernels

using CoordinateTransformations
using Distances
using SpecialFunctions
using StaticArrays

"""
    Kernels.PSFKernel{T} <: AbstractMatrix{T}
"""
abstract type PSFKernel{T} <: AbstractMatrix{T} end

# always inbounds
Base.checkbounds(::Type{Bool}, ::PSFKernel, idx...) = true
Base.checkbounds(::Type{Bool}, ::PSFKernel, idx::CartesianIndex) = true

function indices_from_extent(pos, fwhm, maxsize)
    halfextent = @. maxsize * fwhm / 2
    lower = @. floor(Int, pos - halfextent)
    upper = @. ceil(Int, pos + halfextent)
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

end # module PSFKernels

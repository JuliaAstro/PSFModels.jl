
@doc raw"""
    PSFModels.Gaussian([T=Float64]; fwhm, x=0, y=0, amp=1, maxsize=3, extent=maxsize .* fwhm)
    PSFModels.Gaussian([T=Float64]; fwhm, pos=(0, 0), amp=1, maxsize=3, extent=maxsize .* fwhm)
    PSFModels.Gaussian([T=Float64]; fwhm, r=0, theta=0, origin=(0, 0), amp=1, maxsize=3, extent=maxsize .* fwhm)

An unnormalized bivariate Gaussian distribution. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the model is placed at the origin. The position can also be given as a polar coordinate using `r`/`ρ` and `theta`/`θ`, optionally centered around `origin`.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). For efficient calculations, we recommend using [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis. The `extent` defines the bounding box for the model and is used for the default rendering size.

# Functional form
```
f(x | x̂, FWHM) = exp[-4ln(2) * ||x - x̂|| / FWHM^2]
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a scalar, the Gaussian distribution will be isotropic. If `FWHM` is a vector or tuple, the weighting is applied along each axis (diagonal).
"""
struct Gaussian{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFModel{T}
    fwhm::FT
    pos::VT
    amp::T
    indices::IT

    function Gaussian(fwhm::FT, pos::VT, amp::T, indices::IT) where {T,FT,VT<:AbstractVector,IT<:Tuple}
        new{T,FT,VT,IT}(fwhm, pos, amp, indices)
    end
end

# Alias Normal -> Gaussian
"""
PSFModels.Normal

An alias for [`PSFModels.Gaussian`](@ref)
"""
const Normal = Gaussian

## constructors
# default type is Float64
Gaussian(; kwargs...) = Gaussian(Float64; kwargs...)
function Gaussian(T; fwhm, amp=one(T), maxsize=3, extent = maxsize .* fwhm, position...)
    # get the position from keyword distpatch
    pos = _position(; position...)
    return Gaussian(fwhm, pos, convert(T, amp), indices_from_extent(pos, extent))
end

Base.size(g::Gaussian) = map(length, g.indices)
Base.axes(g::Gaussian) = g.indices

# Gaussian pre-factor for normalizing the exponential
const GAUSS_PRE = -4 * log(2)

# covers isotropic case
function (g::Gaussian{T})(point::AbstractVector) where T
    Δ = sqeuclidean(point, g.pos)
    val = g.amp * exp(GAUSS_PRE * Δ / g.fwhm^2)
    return convert(T, val)
end
# covers vector case
function (g::Gaussian{T,<:Union{Tuple,AbstractVector}})(point::AbstractVector) where T
    weights = SA[1/first(g.fwhm)^2, 1/last(g.fwhm)^2] # manually invert
    Δ = wsqeuclidean(point, g.pos, weights)
    val = g.amp * exp(GAUSS_PRE * Δ)
    return convert(T, val)
end

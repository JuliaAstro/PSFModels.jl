
@doc raw"""
    PSFModels.AiryDisk(fwhm; maxsize=3)
    PSFModels.AiryDisk(position, fwhm; maxsize=3)
    PSFModels.AiryDisk(x, y, fwhm; maxsize=3)
    PSFModels.AiryDisk(::Polar, fwhm; maxsize=3, origin=(0, 0))
    PSFModels.AiryDisk{T}(args...; kwargs...)

An unnormalized Airy disk. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the model is placed at the origin. The position can also be given as a `CoordinateTransformations.Polar`, optionally centered around `origin`.

The `fwhm` can be a scalar (isotropic) or vector/tuple (diagonal). For efficient calculations, we recommend using [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis.

The output type can be specified, and will default to `Float64`. The amplitude is unnormalized, meaning the maximum value will always be 1. This means the models act like a transmission weighting.

# Functional form

The Airy disk is a distribution over the radius `r` (the square-Euclidean distance)

```
f(x | x̂, FWHM) = [ 2J₁(q) / q ]^2
```
where `J₁` is the first order Bessel function of the first kind and
```
q ≈ π * r / (0.973 * FWHM)
```
"""
struct AiryDisk{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFModel{T}
    fwhm::FT
    pos::VT
    amp::T
    indices::IT

    function AiryDisk(::Type{T}, fwhm::FT, pos::VT, amp::T, indices::IT) where {T,VT<:AbstractVector,FT,IT<:Tuple}
        new{T,FT,VT,IT}(fwhm, pos, amp, indices)
    end
end

## constructors
# default type is Float64
AiryDisk(; kwargs...) = AiryDisk(Float64; kwargs...)
function AiryDisk(T; fwhm, amp=one(T), maxsize=3, extent = maxsize .* fwhm, position...)
    # get the position from keyword distpatch
    pos = _position(; position...)
    return AiryDisk(T, fwhm, pos, convert(T, amp), indices_from_extent(pos, extent))
end

Base.size(a::AiryDisk) = map(length, a.indices)
Base.axes(a::AiryDisk) = a.indices

const rz = 3.8317059702075125 / π

function (a::AiryDisk{T})(point::AbstractVector) where T
    radius = a.fwhm * 1.18677
    Δ = euclidean(point, a.pos)
    r = Δ / (radius / rz)
    val = ifelse(iszero(r), a.amp, a.amp * (2 * besselj1(π * r) / (π * r))^2)
    return convert(T, val)
end

function (a::AiryDisk{T,<:Union{AbstractVector,Tuple}})(point::AbstractVector) where T
    weights = SA[(rz / (first(a.fwhm) * 1.18677))^2, (rz / (last(a.fwhm) * 1.18677))^2]
    r = weuclidean(point, a.pos, weights)
    val = ifelse(iszero(r), a.amp, a.amp * (2 * besselj1(π * r) / (π * r))^2)
    return convert(T, val)
end

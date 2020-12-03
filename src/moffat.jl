
@doc raw"""
    PSFModels.Moffat(fwhm; maxsize=3)
    PSFModels.Moffat(position, fwhm; maxsize=3)
    PSFModels.Moffat(x, y, fwhm; maxsize=3)
    PSFModels.Moffat(::Polar, fwhm; maxsize=3, origin=(0, 0))
    PSFModels.Moffat{T}(args...; kwargs...)

An unnormalized Airy disk. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the model is placed at the origin. The position can also be given as a `CoordinateTransformations.Polar`, optionally centered around `origin`.

The `fwhm` can be a scalar (isotropic) or vector/tuple (diagonal). For efficient calculations, we recommend using [StaticArrys](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis.

The output type can be specified, and will default to `Float64`. The amplitude is unnormalized, meaning the maximum value will always be 1. This means the models act like a transmission weighting.

# Functional form
```
f(x | x̂, FWHM) = 1 / [1 + ||x - x̂|| / (FWHM / 2)^2]
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a vector or tuple, the weighting is applied along each axis.
"""
struct Moffat{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFModel{T}
    pos::VT
    fwhm::FT
    indices::IT

    Moffat{T}(pos::VT, fwhm::FT, indices::IT) where {T,VT<:AbstractVector,FT,IT<:Tuple} = new{T,FT,VT,IT}(pos, fwhm, indices)
end

## constructors
# default type is Float64
Moffat(args...; kwargs...) = Moffat{Float64}(args...; kwargs...)
# parse indices from maxsize
Moffat{T}(pos::AbstractVector, fwhm; maxsize=3) where {T} = Moffat{T}(pos, fwhm, indices_from_extent(pos, fwhm, maxsize))
# default position is [0, 0]
Moffat{T}(fwhm; kwargs...) where {T} = Moffat{T}(SA[0, 0], fwhm; kwargs...)
# # parse position to vector
Moffat{T}(x::Number, y::Number, fwhm; kwargs...) where {T} = Moffat{T}(SA[x, y], fwhm; kwargs...)
Moffat{T}(xy::Tuple, fwhm; kwargs...) where {T} = Moffat{T}(SVector(xy), fwhm; kwargs...)
# # translate polar coordinates to cartesian, optionally recentering
Moffat{T}(p::Polar, fwhm; origin=SA[0, 0], kwargs...) where {T} = Moffat(CartesianFromPolar()(p) .+ origin, fwhm; kwargs...)

Base.size(m::Moffat) = map(length, m.indices)
Base.axes(m::Moffat) = m.indices

# scalar case
function (m::Moffat{T})(point::AbstractVector) where T
    hwhm = m.fwhm / 2
    Δ = sqeuclidean(point, m.pos)
    val = inv(1 + Δ / hwhm^2)
    return convert(T, val)
end

# vector case
function (m::Moffat{T,<:Union{AbstractVector,Tuple}})(point::AbstractVector) where T
    weights = SA[(2 / first(m.fwhm))^2, (2 / last(m.fwhm))^2]
    Δ = wsqeuclidean(point, m.pos, weights)
    val = inv(1 + Δ)
    return convert(T, val)
end


@doc raw"""
    PSFModels.Moffat([T=Float64]; fwhm, x=0, y=0, amp=1, alpha=1, maxsize=3, extent=maxsize .* fwhm)
    PSFModels.Moffat([T=Float64]; fwhm, pos, amp=1, alpha=1, maxsize=3, extent=maxsize .* fwhm)
    PSFModels.Moffat([T=Float64]; fwhm, r, theta, amp=1, alpha=1, maxsize=3, extent=maxsize .* fwhm)

Two dimensional Moffat model. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the model is placed at the origin. The position can also be given as a polar coordinate using `r`/`ρ` and `theta`/`θ`, optionally centered around `origin`.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). For efficient calculations, we recommend using [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis. The `extent` defines the bounding box for the model and is used for the default rendering size.

# Functional form
```
f(x | x̂, FWHM, α) = A / (1 + ||x - x̂|| / (FWHM / 2)^2)^α
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a vector or tuple, the weighting is applied along each axis.
"""
struct Moffat{T,FT,VT<:AbstractVector,IT<:Tuple} <: PSFModel{T}
    fwhm::FT
    pos::VT
    amp::T
    alpha::Int
    indices::IT

    function Moffat(::Type{T}, fwhm::FT, pos::VT, amp::T, alpha, indices::IT) where {T,FT,VT<:AbstractVector,IT<:Tuple}
        new{T,FT,VT,IT}(fwhm, pos, amp, alpha, indices)
    end
end

## constructors
# default type is Float64
Moffat(; kwargs...) = Moffat(Float64; kwargs...)
function Moffat(T; fwhm, amp=one(T), alpha=1, maxsize=3, extent = maxsize .* fwhm, position...)
    # get the position from keyword distpatch
    pos = _position(; position...)
    return Moffat(T, fwhm, pos, convert(T, amp), alpha, indices_from_extent(pos, extent))
end


Base.size(m::Moffat) = map(length, m.indices)
Base.axes(m::Moffat) = m.indices

# scalar case
function (m::Moffat{T})(point::AbstractVector) where T
    hwhm = m.fwhm / 2
    Δ = sqeuclidean(point, m.pos)
    val = m.amp / (1 + Δ / hwhm^2)^m.alpha
    return convert(T, val)
end

# vector case
function (m::Moffat{T,<:Union{AbstractVector,Tuple}})(point::AbstractVector) where T
    weights = SA[(2 / first(m.fwhm))^2, (2 / last(m.fwhm))^2]
    Δ = wsqeuclidean(point, m.pos, weights)
    val = m.amp  / (1 + Δ)^m.alpha
    return convert(T, val)
end

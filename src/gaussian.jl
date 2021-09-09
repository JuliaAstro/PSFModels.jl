
@doc raw"""
    PSFModels.Gaussian([T=Float64]; fwhm, x=0, y=0, amp=1, maxsize=3, extent=maxsize .* fwhm)
    PSFModels.Gaussian([T=Float64]; fwhm, pos=(0, 0), amp=1, maxsize=3, extent=maxsize .* fwhm)
    PSFModels.Gaussian([T=Float64]; fwhm, r=0, theta=0, origin=(0, 0), amp=1, maxsize=3, extent=maxsize .* fwhm)

An unnormalized bivariate Gaussian distribution. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the model is placed at the origin. The position can also be given as a polar coordinate using `r`/`ρ` and `theta`/`θ`.

The `fwhm` can be a scalar (isotropic), vector/tuple (diagonal), or a matrix (correlated). For efficient calculations, we recommend using [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis. The `extent` defines the bounding box for the model and is used for the default rendering size.

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

    Gaussian(fwhm::FT, pos::VT, amp::T, indices::IT) where {T,FT,VT<:AbstractVector,IT<:Tuple} = new{T,FT,VT,IT}(fwhm,pos, amp, indices)
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
# parse indices from maxsize
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

@doc raw"""
    PSFModels.Gaussian([T=Float64], X::AbstractVector; kwargs...)

The following canned models are predefined depending on the length of the input vector. These are meant to be used in optimization problems.

| N | Parameters                  |
|---|:----------------------------|
| 1 | fwhm                        |
| 2 | fwhm\_x, fwhm\_y            |
| 3 | x, y, fwhm                  |
| 4 | x, y, fwhm, amp             |
| 5 | x, y, fwhm\_x, fwhm\_y, amp |
"""
Gaussian(T, X::AbstractVector; kwargs...) = Gaussian(T, X, length(X); kwargs...)
Gaussian(X::AbstractVector; kwargs...) = Gaussian(Float64, X; kwargs...)

# option 1: ifelse
function Gaussian(T, X::AbstractVector, N; kwargs...)
    if N == 1
        Gaussian(T; fwhm=X[begin], kwargs...)
    elseif N == 2
        Gaussian(T; fwhm=X, kwargs...)
    elseif N == 3
        Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=X[begin + 2], kwargs...)
    elseif N == 4
        Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=X[begin + 2], amp=X[end], kwargs...)
    elseif N == 5
        Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=@view(X[begin + 2:begin + 3]), amp=X[end], kwargs...)
    end
end

# # option 2: Val lookup
# Gaussian(T, X::AbstractVector, N) = Gaussian(T, X, Val(N))
# Gaussian(T, X::AbstractVector, ::Val{1}) = Gaussian(T; fwhm=X[begin])
# Gaussian(T, X::AbstractVector, ::Val{2}) = Gaussian(T; fwhm=X)
# Gaussian(T, X::AbstractVector, ::Val{3}) = Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=X[begin + 2])
# Gaussian(T, X::AbstractVector, ::Val{4}) = Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=@view(X[begin + 2:end]))
# Gaussian(T, X::AbstractVector, ::Val{5}) = Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=@view(X[begin + 2:begin + 3]), amp=X[end])

# # option 3: dict closure lookup
# const GAUSS_MODEL_DICT = Dict(
#     1 => (T, X) -> Gaussian(T; fwhm=X[begin]),
#     2 => (T, X) -> Gaussian(T; fwhm=X),
#     3 => (T, X) -> Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=X[begin + 2]),
#     4 => (T, X) -> Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=@view(X[begin + 2:end])),
#     5 => (T, X) -> Gaussian(T; x=X[begin], y=X[begin + 1], fwhm=@view(X[begin + 2:begin + 3]), amp=X[end]),
# )
# Gaussian(T, X::AbstractVector, N) = GAUSS_MODEL_DICT[N](T, X)

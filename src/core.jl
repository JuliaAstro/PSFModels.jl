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


# in general, parse to static vector
(model::PSFModel)(point...) = model(SVector(point))
(model::PSFModel)(point::Tuple) = model(SVector(point))
# disallow index to avoid confusion 
(model::PSFModel)(::CartesianIndex) = error("PSF models should be indexed using `getindex`, (equivalently `[]`)")

# getindex just calls model with reversed indices
Base.getindex(model::PSFModel, idx::Vararg{<:Integer,2}) = model(reverse(idx))
# always inbounds
Base.checkbounds(::Type{Bool}, ::PSFModel, idx) = true
Base.checkbounds(::Type{Bool}, ::PSFModel, idx...) = true

# broadcasting hack to slurp other axes (doesn't work for numbers)
Broadcast.combine_axes(psf::PSFModel, other::AbstractMatrix) = axes(other)
Broadcast.combine_axes(other::AbstractMatrix, psf::PSFModel) = axes(other)

## get the position depending on the keyword inputs
function _position(nt::NamedTuple{(:r, :theta)})
    # note: theta is in degrees
    th_rad = deg2rad(nt.theta)
    return CartesianFromPolar()(Polar(nt.r, th_rad))
end
_position(nt::NamedTuple{(:pos,)}) = SVector(nt.pos)
_position(nt::NamedTuple{(:x, :y)}) = SA[nt.x, nt.y]
_position(nt::NamedTuple{()}) = SA[0, 0]

@kwcall _position(r, theta)
@kwcall _position(pos)
@kwcall _position(x=0, y=0)
# theta is in degrees
@kwalias _position [
    ρ => r,
    θ => theta
]

function indices_from_extent(pos, extent)
    halfextent = @. 0.5 * extent
    lower = @. round(Int, pos - halfextent)
    upper = @. round(Int, pos + halfextent)
    return last(lower):last(upper), first(lower):first(upper)
end
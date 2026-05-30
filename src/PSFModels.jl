module PSFModels

import ADTypes
import ForwardDiff
import Optim
using ConstructionBase: constructorof, getfields, setproperties # Use these to query / update structs in a generic way for fitting
using Rotations: RotMatrix
using SpecialFunctions: besselj1
using StaticArrays: SA

export gaussian, normal, airydisk, moffat
export GaussianPSFSymmetric, evaluate, centroid, integral, fit_deriv, fit_hessian

const BivariateLike = Union{<:Tuple,<:AbstractVector}

function rotate_point(dx, dy, theta)
    # generate rotation matrix
    # (theta is degrees CCW from x-axis)
    R = RotMatrix{2}(-deg2rad(theta))
    # rotate points
    return R * SA[dx, dy]
end

"""Abstract type for PSF models. All PSF models should be subtypes of this abstract type, and implement the following methods:"""
abstract type AbstractPSFModel{T} end
Base.Broadcast.broadcastable(m::AbstractPSFModel) = Ref(m)

"""
    evaluate(model::AbstractPSFModel{T}, x::Real, y::Real)::T
Evaluate the PSF model at position `(x, y)`.
"""
function evaluate end

"""
    centroid(model::AbstractPSFModel{T})::Tuple{T, T}
Return the centroid of the PSF model as a tuple `(x, y)`;
default implementation assumes the centroid is given by fields `x` and `y` in the model struct.
"""
function centroid(model::AbstractPSFModel)
    if hasproperty(model, :x) && hasproperty(model, :y)
        return (model.x, model.y)
    else
        error("Model does not have `x` and `y` fields; either add them or implement `centroid(model)` for this model type.")
    end
end

"""
    integral(model::AbstractPSFModel{T})::T
Return the integral of the PSF model over all space;
default implementation assumes the integral is given by a field `flux` in the model struct.
"""
integral(model::AbstractPSFModel) = hasproperty(model, :flux) ? model.flux : error("Model does not have a `flux` field; either add one or implement `integral(model)` for this model type.")

"""
    fit_deriv(model::AbstractPSFModel{T}, x::Real, y::Real)
Returns the partial derivatives of the `model`
with respect to the parameters.
"""
function fit_deriv end

"""
    fit_hessian(model::AbstractPSFModel{T}, x::Real, y::Real)
Returns the second partial derivatives of the `model`
with respect to the parameters.
"""
function fit_hessian end

include("gaussian.jl")
include("airy.jl")
include("moffat.jl")
include("plotting.jl")

# codegen for common functionality
# if you add a new model, make sure it gets added here
for model in (:gaussian, :airydisk, :moffat)
    @eval begin
        $model(px, py; kwargs...) = $model(Float64, px, py; kwargs...)
        $model(point::BivariateLike; kwargs...) = $model(point...; kwargs...)
        $model(T, point::BivariateLike; kwargs...) = $model(T, point...; kwargs...)
        $model(idx::CartesianIndex; kwargs...) = $model(idx.I; kwargs...)
        $model(T, idx::CartesianIndex; kwargs...) = $model(T, idx.I; kwargs...)
        $model(; kwargs...) = (point...) -> $model(_curried_point(point...); kwargs...)
        $model(::Type{T}; kwargs...) where {T} = (point...) -> $model(T, _curried_point(point...); kwargs...)
    end
end

_curried_point(P::BivariateLike) = P
_curried_point(point...) = Tuple(point)

include("fitting.jl")

end # module PSFModels

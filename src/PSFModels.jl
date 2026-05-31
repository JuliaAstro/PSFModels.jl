module PSFModels

import ADTypes
import ForwardDiff
using LinearAlgebra: inv
import LossFunctions
import NLSolversBase
import Optim
using ConstructionBase: constructorof, getfields, getproperties, setproperties # Use these to query / update structs in a generic way for fitting
using Rotations: RotMatrix
using SpecialFunctions: besselj1
using StaticArrays: SA

export gaussian, normal, airydisk, moffat
export CircularGaussianPSF, GaussianPSF, evaluate, evaluate_fg, evaluate_fgh, centroid, integral

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
(model::AbstractPSFModel)(x, y) = evaluate(model, x, y)
(model::AbstractPSFModel)(idx::CartesianIndex) = evaluate(model, Tuple(idx)...)
evaluate(model::AbstractPSFModel, idx::CartesianIndex) = evaluate(model, Tuple(idx)...)

"""
    evaluate(model::AbstractPSFModel{T}, x::Real, y::Real)::T
Evaluate the PSF model at position `(x, y)`.
"""
function evaluate end

"""
    centroid(model::AbstractPSFModel{T}) → (x::T, y::T)
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
    fwhm(model::AbstractPSFModel) → (x_fwhm::T, y_fwhm::T)
Return the full width at half maximum (FWHM) of the PSF model as a tuple `(x_fwhm, y_fwhm)` in the x and y directions. By default, this function checks for a single `fwhm` field and returns it for both axes, or separate `x_fwhm` and `y_fwhm` fields if they exist. Models with more complex definitions of FWHM should implement their own version of this function.
"""
function fwhm(model::AbstractPSFModel)
    if hasproperty(model, :fwhm)
        return (model.fwhm, model.fwhm)
    elseif hasproperty(model, :x_fwhm) && hasproperty(model, :y_fwhm)
        return (model.x_fwhm, model.y_fwhm)
    else
        error("Model does not have `fwhm` or `x_fwhm` and `y_fwhm` fields; either add them or implement `fwhm(model)` for this model type.")
    end
end

"""
    evaluate_fg(model::AbstractPSFModel{T}, x::Real, y::Real) → (f::T, G::SVector{T})
Returns the model value `f` and partial derivatives of the `model`
with respect to the parameters `G` at position `(x, y)`.
"""
function evaluate_fg end

"""
    evaluate_fgh(model::AbstractPSFModel{T}, x::Real, y::Real) → (f::T, G::SVector{T}, H::SMatrix{T})
Returns the model value `f`, partial derivatives `G`, and Hessian matrix `H` of the `model`
with respect to the parameters at position `(x, y)`.
"""
function evaluate_fgh end

"""
    ellipse_bounds(a, b, θ)
Returns the bounding box of an ellipse with semi-major axis `a`, 
semi-minor axis `b`, and rotation angle `θ` (degrees CCW from x-axis).

```jldoctest
julia> using PSFModels: ellipse_bounds

julia> ellipse_bounds(3, 2, 0)
(3.0, 2.0)

julia> ellipse_bounds(3, 2, 90)
(2.0, 3.0)

julia> round.(ellipse_bounds(3, 2, 45); digits=3)
(2.55, 2.55)
```
"""
function ellipse_bounds(a, b, θ)
    θ = deg2rad(θ)
    cosθ = cos(θ)
    sinθ = sin(θ)
    x_bound = hypot(a * cosθ, b * sinθ)
    y_bound = hypot(a * sinθ, b * cosθ)
    return x_bound, y_bound
end

function bounding_box(model::AbstractPSFModel, factor=5)
    # default bounding box is 5x5 around centroid, but specific models can override this
    x0, y0 = centroid(model)
    FWHM = fwhm(model)
    a, b = factor * FWHM[1], factor * FWHM[2]
    dx, dy = ellipse_bounds(a, b, 0)
    return (x0 - 2.5, x0 + 2.5), (y0 - 2.5, y0 + 2.5)
end

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
include("fitting_struct.jl")

end # module PSFModels

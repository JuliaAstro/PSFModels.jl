module PSFModels

import ADTypes
import ForwardDiff
using LinearAlgebra: inv, cholesky!, ldiv!, I, Symmetric, pinv
import LossFunctions
import NLSolversBase
import Optim
using ConstructionBase: constructorof, getfields, getproperties, setproperties # Use these to query / update structs in a generic way for fitting
using Rotations: RotMatrix
using SpecialFunctions: besselj, besselj0, besselj1, erf
using StaticArrays: SA, SVector, MVector, MMatrix
using Statistics: median, mean

export gaussian, normal, airydisk, moffat
export AiryPSF, CircularGaussianPSF, GaussianPSF, CircularGaussianPRF, GaussianPRF, CircularMoffat, MoffatPSF, evaluate, centroid, integral, render, render!, peak, amplitude, effective_area

"""Abstract type for PSF models. All PSF models should be subtypes of this abstract type, and implement the following methods:"""
abstract type AbstractPSFModel{T} end
Base.Broadcast.broadcastable(m::AbstractPSFModel) = Ref(m)
(model::AbstractPSFModel)(x, y) = evaluate(model, x, y)
(model::AbstractPSFModel)(idx::CartesianIndex) = evaluate(model, Tuple(idx)...)
evaluate(model::AbstractPSFModel, idx::CartesianIndex) = evaluate(model, Tuple(idx)...)
evaluate_fg(model::AbstractPSFModel, idx::CartesianIndex) = evaluate_fg(model, Tuple(idx)...)
evaluate_fgh(model::AbstractPSFModel, idx::CartesianIndex) = evaluate_fgh(model, Tuple(idx)...)
function evaluate_fg(model::AbstractPSFModel, x, y, free_idx::AbstractVector)
    f, g = evaluate_fg(model, x, y)
    return f, view(g, free_idx)
end
function evaluate_fg(model::AbstractPSFModel, x, y, free_idx::SVector)
    f, g = evaluate_fg(model, x, y)
    return f, g[free_idx]
end

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
function integral(model::AbstractPSFModel)
    if hasproperty(model, :flux)
        return model.flux
    else
        error("Model does not have a `flux` field; either add one or implement `integral(model)` for this model type.")
    end
end

"""
    background(model::AbstractPSFModel{T})::T
Return the background level of the PSF model; 
if `bkg` field exists, return that, otherwise return `zero(T)`.
"""
function background(model::AbstractPSFModel)
    if hasproperty(model, :bkg)
        return model.bkg
    else
        error("Model does not have a `bkg` field; either add one or implement `background(model)` for this model type.")
    end
end

"""
    peak(model::AbstractPSFModel{T})::T

Return the peak value of the PSF model. By default, this
function evaluates the model at its centroid, but
models can override this.
"""
peak(model::AbstractPSFModel) = evaluate(model, centroid(model)...)

"""
    amplitude(model::AbstractPSFModel{T})::T

Return the amplitude of the PSF model, defined as the peak value
minus the background. Default implementation is
`peak(model) - background(model)`.
"""
amplitude(model::AbstractPSFModel) = peak(model) - background(model)

"""
    effective_area(model::AbstractPSFModel{T})::T

Return the effective area of the PSF model, defined as 
``\\frac{\\left(\\int PSF(x, y) dx dy\\right)^2}{\\int PSF(x, y)^2 dx dy}``. 
For PSF fitting
photometry, this is the effective number of noisy pixels that contribute
to the measurement -- the variance of the flux measurement is approximately
the variance of the background noise per pixel times the effective area. 
Models with more complex definitions of effective area should implement 
their own version of this function.

For a star image with PSF `P` and background noise per pixel with variance `σ²`, 
the maximum likelihood estimate of the flux is ``F^\\hat = \\frac{\\sum_i P_i (D_i - B_i)}{\\sum_i P_i^2}``, where `D_i` is the observed data, `B_i` is the background, and `P_i` is the PSF value at pixel `i`. The variance of the flux measurement is ``var(F^\\hat) = σ² * effective_area(P)``.
"""
function effective_area(model::AbstractPSFModel) end

"""
    fwhm(model::AbstractPSFModel{T}) → (x_fwhm::T, y_fwhm::T)

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
    theta(model::AbstractPSFModel{T}) → θ::T

Return the rotation angle `theta` of the PSF model in degrees CCW from the x-axis. By default, this function checks for a `theta` field and returns it, or **returns zero if no such field exists**. Models with more complex definitions of rotation should implement their own version of this function.
"""
function theta(model::AbstractPSFModel{T}) where {T}
    if hasproperty(model, :theta)
        return model.theta
    else
        return zero(T)
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
    ellipse_bounds(a, b, θ) → (x_bound, y_bound)

Returns the half-width and half-height of the smallest axis-aligned
rectangle enclosing an ellipse with semi-major axis `a`, semi-minor
axis `b`, and rotation angle `θ` (degrees counterclockwise from the
x-axis).

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
    sinθ, cosθ = sincos(θ)
    x_bound = hypot(a * cosθ, b * sinθ)
    y_bound = hypot(a * sinθ, b * cosθ)
    return x_bound, y_bound
end

"""
    extent(model::AbstractPSFModel, fwhm_factor=5) → (x_range::Tuple, y_range::Tuple)

Returns the extent of the PSF model which is typically useful for fitting, plotting, etc., 
`((x_min, x_max), (y_min, y_max))`. By default, the extent is the smallest axis-aligned
rectangle enclosing an ellipse centered on the model centroid whose major and minor axis
lengths are `fwhm_factor` times the model FWHM values. Models with more
complex shapes can override this function to provide a more appropriate extent.

```jldoctest
julia> using PSFModels: extent, CircularGaussianPSF, GaussianPSF

julia> extent(CircularGaussianPSF(x=10, y=20, fwhm=5, flux=30, bkg=1), 5)
((-2.5, 22.5), (7.5, 32.5))

julia> extent(GaussianPSF(x=10, y=20, x_fwhm=5, y_fwhm=3, theta=90, flux=30, bkg=1), 5)
((2.5, 17.5), (7.5, 32.5))
```
"""
function extent(model::AbstractPSFModel, fwhm_factor = 5)
    # default extent is 5x5 around centroid, but specific models can override this
    x0, y0 = centroid(model)
    FWHM = fwhm(model)
    a, b = fwhm_factor * FWHM[1] / 2, fwhm_factor * FWHM[2] / 2
    dx, dy = ellipse_bounds(a, b, theta(model))
    return (x0 - dx, x0 + dx), (y0 - dy, y0 + dy)
end

"""
    render!(out::AbstractMatrix, model::AbstractPSFModel, inds)

Fill the pre-allocated matrix `out` with `evaluate(model, px, py)` for each pixel
`(px, py)` in `inds`. The first axis of `out` corresponds to `inds[1]` and the
second axis to `inds[2]`, so offset ranges are handled correctly.
"""
function render!(out::AbstractMatrix, model::AbstractPSFModel, inds)
    xs, ys = inds
    @inbounds for i in eachindex(xs)
        px = xs[i]
        for j in eachindex(ys)
            out[i, j] = evaluate(model, px, ys[j])
        end
    end
    return out
end

"""
    render(model::AbstractPSFModel{T}, inds) → Matrix{float(T)}

Allocate and return a matrix whose `[i, j]` entry is `evaluate(model, px, py)` for
`px = inds[1][i]` and `py = inds[2][j]`.
"""
function render(model::AbstractPSFModel{T}, inds) where {T}
    out = Matrix{float(T)}(undef, length(inds[1]), length(inds[2]))
    return render!(out, model, inds)
end

"""
    render(model::AbstractPSFModel) → Matrix

Allocate and return a matrix covering the region returned by `extent(model)`,
rounding the bounds to the nearest integer.
"""
function render(model::AbstractPSFModel)
    (x_lo, x_hi), (y_lo, y_hi) = extent(model)
    inds = (round(Int, x_lo):round(Int, x_hi), round(Int, y_lo):round(Int, y_hi))
    return render(model, inds)
end

include("parametric_models.jl")
include("functional_models.jl")
include("plotting.jl")

######################################################
# Fitting code

include("fitting.jl")
include("fitting_struct.jl")
include("levenberg_marquardt.jl")

end # module PSFModels

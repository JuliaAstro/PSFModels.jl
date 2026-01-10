
@doc raw"""
    gaussian([T=Float64], point; x, y, fwhm, amp=1, theta=0, bkg=0)
    gaussian([T=Float64], px, py; x, y, fwhm, amp=1, theta=0, bkg=0)

An unnormalized bivariate Gaussian distribution. The position can be specified
in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate
arguments. If `theta` is given, the PSF will be rotated by `theta` degrees
counter-clockwise from the x-axis. If `bkg` is given it will be added as a
scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in
mind that `theta` has no effect for isotropic distributions and is degenerate
with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm`
tuple)

# Functional form
```math
f(x | x̂, \mathrm{FWHM}) = \exp[-4 \ln(2) ⋅ ||x - x̂|| / \mathrm{FWHM}^2]
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the
square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a
scalar, the Gaussian distribution will be isotropic. If `FWHM` is a vector or
tuple, the weighting is applied along each axis (diagonal).
"""
gaussian(T, px, py; x, y, fwhm, amp=one(T), theta=0, bkg=0) =
    convert(T, _gaussian(px, py, x, y, fwhm, amp, theta, bkg))


"""
    normal

An alias for [`gaussian`](@ref)
"""
const normal = gaussian

# this is the factor to convert 1/(2σ²) to 1/(2fwhm²)
const GAUSS_PRE = -4 * log(2)

# isotropic
function _gaussian(px, py, x, y, fwhm, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic gaussian is not affected by non-zero rotation angle $theta"
    # unnormalized gaussian likelihood
    sqmahab = dx^2 + dy^2
    sqmahab /= fwhm^2
    return amp * exp(GAUSS_PRE * sqmahab) + background
end

# bivariate
function _gaussian(px, py, x, y, fwhm::BivariateLike, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    if !iszero(theta)
        dx, dy = rotate_point(dx, dy, theta)
    end
    # unnormalized gaussian likelihood
    fwhmx, fwhmy = fwhm
    sqmahab = (dx / fwhmx)^2 + (dy / fwhmy)^2
    return amp * exp(GAUSS_PRE * sqmahab) + background
end

## gradients

# isotropic
function fgrad(g::typeof(gaussian), point::AbstractVector)
    f = g(point)

    xdiff = first(point) - first(g.pos)
    ydiff = last(point) - last(g.pos)
    dfdpos = -2 * GAUSS_PRE * f / g.fwhm^2 .* SA[xdiff, ydiff]
    dfdfwhm = -2 * GAUSS_PRE * f * (xdiff^2 + ydiff^2) / g.fwhm^3
    dfdamp = f / g.amp
    return f, dfdpos, dfdfwhm, dfdamp
end

# diagonal
#function fgrad(g::Gaussian{T,<:Union{Tuple,AbstractVector}}, point::AbstractVector) where T
#    f = g(point)
#
#    xdiff = first(point) - first(g.pos)
#    ydiff = last(point) - last(g.pos)
#    dfdpos = -2 * GAUSS_PRE * f .* SA[xdiff / first(g.fwhm)^2, ydiff / last(g.fwhm)^2]
#    dfdfwhm = -2 * GAUSS_PRE * f .* SA[xdiff^2 / first(g.fwhm)^3, ydiff^2 / last(g.fwhm)^3]
#    dfda = f / g.amp
#    return f, dfdpos, dfdfwhm, dfda
#end

function frule((Δpsf, Δp), g::typeof(gaussian), point::AbstractVector)
    f, dfdpos, dfdfwhm, dfda = fgrad(g, point)
    Δf = dot(dfdpos, Δpsf.pos) + dot(dfdfwhm, Δpsf.fwhm) + dfda * Δpsf.amp
    Δf -= dot(dfdpos, Δp)
    return f, Δf
end

function rrule(g::typeof(gaussian), point::AbstractVector)
    f, dfdpos, dfdfwhm, dfda = fgrad(g, point)
    function Gaussian_pullback(Δf)
        ∂pos = dfdpos .* Δf
        ∂fwhm = dfdfwhm .* Δf
        ∂g = Tangent{G}(pos=∂pos, fwhm=∂fwhm, amp=dfda * Δf, indices=NoTangent())
        ∂pos = dfdpos .* -Δf
        return ∂g, ∂pos
    end
    return f, Gaussian_pullback
end

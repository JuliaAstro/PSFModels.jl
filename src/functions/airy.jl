
@doc raw"""
    PSFModels.airydisk([T=Float64], point; x, y, fwhm, amp=1, theta=0)
    PSFModels.airydisk([T=Float64], px, py; x, y, fwhm, amp=1, theta=0)

An unnormalized Airy disk. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta` is given, the PSF will be rotated by `theta` degrees counter-clockwise from the x-axis.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in mind that `theta` has no effect for isotropic distributions and is degenerate with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm` tuple)

# Functional form

The Airy disk is a distribution over the radius `r` (the square-Euclidean distance)

```
f(x | x̂, FWHM) = [ 2J₁(q) / q ]^2
```
where `J₁` is the first-order Bessel function of the first kind and
```
q ≈ π * r / (0.973 * FWHM)
```
"""
airydisk(T, px, py; x, y, fwhm, amp=one(T), theta=0) = convert(T, _airydisk(px, py, x, y, fwhm, amp, theta))
airydisk(px, py; kwargs...) = airydisk(Float64, px, py; kwargs...)
airydisk(point::BivariateLike; kwargs...) = airydisk(point...; kwargs...)
airydisk(idx::CartesianIndex; kwargs...) = airydisk(idx.I; kwargs...)

# factor for scaling radius in terms of the fwhm
const rz = 1.18677 * π / 3.8317059702075125

# isotropic
function _airydisk(px, py, x, y, fwhm, amp, theta)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic airydisk is not affected by non-zero rotation angle $theta"
    # unnormalized airydisk
    r = sqrt(dx^2 + dy^2)
    # short-circuit
    iszero(r) && return amp
    r /= (fwhm * rz)
    return amp * (2 * besselj1(π * r) / (π * r))^2
end

# bivariate
function _airydisk(px, py, x, y, fwhm::BivariateLike, amp, theta)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    if !iszero(theta)
        dx, dy = rotate_point(dx, dy, theta)
    end
    # unnormalized airydisk
    fwhmx, fwhmy = fwhm
    r = sqrt((dx / (fwhmx * rz))^2 + (dy / (fwhmy * rz))^2)
    # short-circuit
    iszero(r) && return amp
    return amp * (2 * besselj1(π * r) / (π * r))^2
end

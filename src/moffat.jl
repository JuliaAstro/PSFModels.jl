
@doc raw"""
    moffat([T=Float64], point; x, y, fwhm, alpha=1, amp=1, theta=0, bkg=0)
    moffat([T=Float64], px, py; x, y, fwhm, alpha=1, amp=1, theta=0, bkg=0)

Two dimensional Moffat model. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta` is given, the PSF will be rotated by `theta` degrees counter-clockwise from the x-axis. If `bkg` is given it will be added as a scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in mind that `theta` has no effect for isotropic distributions and is degenerate with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm` tuple)

# Functional form
```
f(x | x̂, FWHM, α) = A / (1 + ||x - x̂|| / (FWHM / 2)^2)^α
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a vector or tuple, the weighting is applied along each axis.
"""
moffat(T, px, py; x, y, fwhm, alpha=1, amp=one(T), theta=0, bkg=0) = convert(T, _moffat(px, py, x, y, fwhm, alpha, amp, theta, bkg))

# isotropic
function _moffat(px, py, x, y, fwhm, alpha, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic moffat is not affected by non-zero rotation angle $theta"
    # unnormalized moffat
    dist = dx^2 + dy^2
    dist /= (fwhm / 2)^2
    return amp / (1 + dist)^alpha + background
end

# bivariate
function _moffat(px, py, x, y, fwhm::BivariateLike, alpha, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    if !iszero(theta)
        dx, dy = rotate_point(dx, dy, theta)
    end
    # unnormalized moffat
    fwhmx, fwhmy = fwhm
    # unnormalized moffat
    dist = (dx / (fwhmx/2))^2 + (dy / (fwhmy/2))^2
    return amp / (1 + dist)^alpha + background
end

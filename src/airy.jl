
@doc raw"""
    airydisk([T=Float64], point; x, y, fwhm, ratio=0, amp=1, theta=0, bkg=0)
    airydisk([T=Float64], px, py; x, y, fwhm, ratio=0, amp=1, theta=0, bkg=0)

An unnormalized Airy disk. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta` is given, the PSF will be rotated by `theta` degrees counter-clockwise from the x-axis. If `bkg` is given it will be added as a scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in mind that `theta` has no effect for isotropic distributions and is degenerate with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm` tuple)

If `ratio` is supplied, this will be the Airy pattern for a centrally-obscured aperture (e.g., a Newtonian telescope). This has a slightly expanded functional form, and in general the central Airy disk will be smaller and the first Airy ring will be brighter.

# Functional form

The Airy disk is a distribution over the radius `r` (the square-Euclidean distance)

```
f(x | x̂, FWHM) = [ 2J₁(q) / q ]^2
```
where `J₁` is the first-order Bessel function of the first kind and
```
q ≈ π * r * D/ λ ≈ π * r / (0.973 * FWHM)
```

If user a non-zero central obscuration via `ratio`, the functional form becomes

```
f(x | x̂, FWHM, ϵ) = [ 2J₁(q) / q - 2ϵJ₁(ϵq) / q ]^2 / (1 - ϵ^2)^2
```
where `ϵ` is the ratio (0 ≤ `ϵ` < 1).
"""
airydisk(T, px, py; x, y, fwhm, ratio=0, amp=one(T), theta=0, bkg=0) = convert(T, _airydisk(px, py, x, y, fwhm, ratio, amp, theta, bkg))

# factor for scaling radius in terms of the fwhm
const AIRY_PRE = π / 0.973

# isotropic
function _airydisk(px, py, x, y, fwhm, ratio, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic airydisk is not affected by non-zero rotation angle $theta"
    # unnormalized airydisk
    r = sqrt(dx^2 + dy^2)
    # short-circuit
    iszero(r) && return amp / (1 - ratio^2)^2 + background
    q = AIRY_PRE * r / fwhm
    I1 = 2 * besselj1(q) / q
    if iszero(ratio)
        return amp * I1^2 + background
    end
    I2 = 2 * ratio * besselj1(q * ratio) / q
    return amp * ((I1 - I2) / (1 - ratio^2))^2 + background
end

# bivariate
function _airydisk(px, py, x, y, fwhm::BivariateLike, ratio, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    if !iszero(theta)
        dx, dy = rotate_point(dx, dy, theta)
    end
    # unnormalized airy disk
    fwhmx, fwhmy = fwhm
    q = AIRY_PRE * sqrt((dx / fwhmx)^2 + (dy / fwhmy)^2)
    # short-circuit
    iszero(q) && return amp / (1 - ratio^2)^2 + background
    I1 = 2 * besselj1(q) / q
    if iszero(ratio)
        return amp * I1^2 + background
    end
    I2 = 2 * ratio * besselj1(q * ratio) / q
    return amp * ((I1 - I2) / (1 - ratio^2))^2 + background
end

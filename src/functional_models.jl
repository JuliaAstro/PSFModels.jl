# Functional model interface

const BivariateLike = Union{<:Tuple,<:AbstractVector}

function rotate_point(dx, dy, theta)
    # generate rotation matrix
    # (theta is degrees CCW from x-axis)
    R = RotMatrix{2}(-deg2rad(theta))
    # rotate points
    return R * SA[dx, dy]
end

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
function gaussian(T, px, py; x, y, fwhm, amp=1, theta=0, bkg=0)
    flux = amp * (π * (fwhm isa BivariateLike ? prod(fwhm) : fwhm^2) / -GAUSS_PRE)
    model = if fwhm isa BivariateLike
        GaussianPSF(x, y, fwhm..., theta, flux, bkg)
    else
        !iszero(theta) && @warn "isotropic gaussian is not affected by non-zero rotation angle $theta"
        CircularGaussianPSF(x, y, fwhm, flux, bkg)
    end
    return convert(T, evaluate(model, px, py))
end

"""
    normal

An alias for [`gaussian`](@ref)
"""
const normal = gaussian


@doc raw"""
    airydisk([T=Float64], point; x, y, fwhm, ratio=0, amp=1, theta=0, bkg=0)
    airydisk([T=Float64], px, py; x, y, fwhm, ratio=0, amp=1, theta=0, bkg=0)

An unnormalized Airy disk. The position can be specified in `(x, y)` coordinates
as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta` is given,
the PSF will be rotated by `theta` degrees counter-clockwise from the x-axis. If
`bkg` is given it will be added as a scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in
mind that `theta` has no effect for isotropic distributions and is degenerate
with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm`
tuple)

If `ratio` is supplied, this will be the Airy pattern for a centrally-obscured
aperture (e.g., a Newtonian telescope). This has a slightly expanded functional
form, and in general the central Airy disk will be smaller and the first Airy
ring will be brighter.

# Functional form

The Airy disk is a distribution over the radius `r` (the square-Euclidean distance)

```math
f(x | x̂, \mathrm{FWHM}) = [ 2J₁(q) / q ]^2
```
where `J₁` is the first-order Bessel function of the first kind and
```math
q ≈ π r D / λ ≈ π r / (0.973 × \mathrm{FWHM})
```

If user a non-zero central obscuration via `ratio`, the functional form becomes

```math
f(x | x̂, \mathrm{FWHM}, ϵ) = [ 2J₁(q) / q - 2ϵJ₁(ϵq) / q ]^2 / (1 - ϵ^2)^2
```
where ``ϵ`` is the ratio (``0 ≤ ϵ < 1``).
"""
airydisk(T, px, py; x, y, fwhm, ratio=0, amp=one(T), theta=0, bkg=0) =
    convert(T, _airydisk(px, py, x, y, fwhm, ratio, amp, theta, bkg))

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
    iszero(r) && return amp + background
    q = AIRY_PRE * r / fwhm
    I1 = 2 * besselj1(q) / q
    iszero(ratio) && return amp * I1^2 + background
    I2 = 2 * ratio * besselj1(q * ratio) / q
    return amp * (I1 - I2)^2 + background
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
    iszero(q) && return amp + background
    I1 = 2 * besselj1(q) / q
    iszero(ratio) && return amp * I1^2 + background
    I2 = 2 * ratio * besselj1(q * ratio) / q
    return amp * (I1 - I2)^2+ background
end


@doc raw"""
    moffat([T=Float64], point; x, y, fwhm, alpha=1, amp=1, theta=0, bkg=0)
    moffat([T=Float64], px, py; x, y, fwhm, alpha=1, amp=1, theta=0, bkg=0)

Two dimensional Moffat model. The position can be specified in `(x, y)`
coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta`
is given, the PSF will be rotated by `theta` degrees counter-clockwise from the
x-axis. If `bkg` is given it will be added as a scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in
mind that `theta` has no effect for isotropic distributions and is degenerate
with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm`
tuple)

# Functional form
```math
f(x | x̂, \mathrm{FWHM}, α) = A (1 + (||x - x̂|| / \mathrm{FWHM} / 2)^2)^{-α}
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the
distance, and `FWHM` is the full width at half-maximum.
If `fwhm` is a vector or tuple, the weighting is applied along each axis.

Note that this function technically uses the half width at half-maximum, defined
as ``\mathrm{HWHM} = \mathrm{FWHM}/2``, but for compatibility with the other
models, `fwhm` is used as an input parameter instead.
"""
moffat(T, px, py; x, y, fwhm, alpha=1, amp=one(T), theta=0, bkg=0) =
    convert(T, _moffat(px, py, x, y, fwhm, alpha, amp, theta, bkg))

# isotropic
function _moffat(px, py, x, y, fwhm, alpha, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic moffat is not affected by non-zero rotation angle $theta"
    # unnormalized moffat
    gamma = _moffat_fwhm_to_gamma(fwhm, alpha)
    dist = (dx / gamma)^2 + (dy / gamma)^2
    return amp / (1 + dist)^alpha + background
end

# http://openafox.com/science/peak-function-derivations.html#moffat
_moffat_gamma_to_fwhm(gamma, alpha) = 2 * gamma * sqrt(2^(1/alpha) - 1)
_moffat_fwhm_to_gamma(fwhm, alpha) = fwhm / (2 * sqrt(2^(1/alpha) - 1))

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
    gammax, gammay = _moffat_fwhm_to_gamma.(fwhm, alpha)
    # unnormalized moffat
    dist = (dx / gammax)^2 + (dy / gammay)^2
    return amp / (1 + dist)^alpha + background
end


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
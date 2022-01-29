
@doc raw"""
    PSFModels.gaussian([T=Float64], point; x, y, fwhm, amp=1, theta=0)
    PSFModels.gaussian([T=Float64], px, py; x, y, fwhm, amp=1, theta=0)

An unnormalized bivariate Gaussian distribution. The position can be specified in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. By default the model is placed at the origin. The position can also be given as a polar coordinate using `r`/`ρ` and `theta`/`θ`, optionally centered around `origin`.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). For efficient calculations, we recommend using [StaticArrays](https://github.com/JuliaArrays/StaticArrays.jl). Here, `maxsize` is a multiple of the fwhm, and can be given as a scalar or as a tuple for each axis. The `extent` defines the bounding box for the model and is used for the default rendering size.

# Functional form
```
f(x | x̂, FWHM) = exp[-4ln(2) * ||x - x̂|| / FWHM^2]
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a scalar, the Gaussian distribution will be isotropic. If `FWHM` is a vector or tuple, the weighting is applied along each axis (diagonal).
"""
gaussian(T, px, py; x, y, fwhm, amp=one(T), theta=0) = convert(T, _gaussian(px, py, x, y, fwhm, amp, theta))
gaussian(px, py; kwargs...) = gaussian(Float64, px, py; kwargs...)
gaussian(point::BivariateLike; kwargs...) = gaussian(point...; kwargs...)
gaussian(idx::CartesianIndex; kwargs...) = gaussian(idx.I; kwargs...)

"""
    PSFModels.normal

An alias for [`PSFModels.gaussian`](@ref)
"""
const normal = gaussian

# this is the factor to convert 1/(2σ²) to 1/(2fwhm²)
const GAUSS_PRE = -4 * log(2)

# isotropic
function _gaussian(px, py, x, y, fwhm, amp, theta)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate 
    !iszero(theta) && @warn "isotropic gaussian is not affected by non-zero rotation angle $theta"
    # unnormalized gaussian likelihood
    sqmahab = dx^2 + dy^2
    sqmahab /= fwhm^2
    return amp * exp(GAUSS_PRE * sqmahab)
end

# bivariate
function _gaussian(px, py, x, y, fwhm::BivariateLike, amp, theta)
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
    return amp * exp(GAUSS_PRE * sqmahab)
end

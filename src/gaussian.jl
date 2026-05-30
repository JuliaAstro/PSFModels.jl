
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

"""
    CircularGaussianPSF(x, y, fwhm, flux, bkg) -> CircularGaussianPSF{T}

Circular, symmetric Gaussian PSF with centroid `(x, y)` and FWHM given by `fwhm`. The `flux` is the integral of the PSF over all space, and `bkg` is a scalar background level added to the PSF.

```jldoctest
julia> using PSFModels: CircularGaussianPSF

julia> CircularGaussianPSF(1.0, 2.0, 3.0, 4.0, 5.0) isa CircularGaussianPSF{Float64}
true
```
"""
Base.@kwdef struct CircularGaussianPSF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    fwhm::T
    flux::T
    bkg::T
    function CircularGaussianPSF(x, y, fwhm, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(fwhm), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T # promote to Float if all inputs are integers
        return new{T}(T(x), T(y), T(fwhm), T(flux), T(bkg))
    end
end

function evaluate(model::CircularGaussianPSF{T}, px, py) where T
    dx = px - model.x
    dy = py - model.y
    sqmahab = (dx^2 + dy^2) / model.fwhm^2
    norm = -(π * model.fwhm^2 / T(GAUSS_PRE))
    amp = model.flux / norm
    return muladd(amp, exp(T(GAUSS_PRE) * sqmahab), model.bkg)
end

function fit_deriv(model::CircularGaussianPSF{T}, px, py) where T
    _gauss_pre = T(GAUSS_PRE)
    dx = px - model.x
    dy = py - model.y
    fwhm = model.fwhm
    sqmahab = (dx^2 + dy^2) / fwhm^2
    norm = -(T(π) * fwhm^2 / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)
    dg_dx = -2 * amp * _gauss_pre * g * dx / fwhm^2
    dg_dy = -2 * amp * _gauss_pre * g * dy / fwhm^2
    dg_dfwhm = amp * g * (-2 / fwhm - 2 * _gauss_pre * sqmahab / fwhm)
    dg_dflux = g / norm
    dg_dbkg = one(T)
    return SA[dg_dx, dg_dy, dg_dfwhm, dg_dflux, dg_dbkg]
end

function fit_hessian(model::CircularGaussianPSF{T}, px, py) where T
    _gauss_pre = T(GAUSS_PRE)
    dx = px - model.x
    dy = py - model.y
    fwhm = model.fwhm
    fwhm² = fwhm^2
    sqmahab = (dx^2 + dy^2) / fwhm²
    norm = -(T(π) * fwhm² / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)

    # Shared sub-expressions
    Ag = amp * g
    γ_f2 = _gauss_pre / fwhm² # γ/fwhm²
    γr2_f = _gauss_pre * sqmahab / fwhm # γr²/fwhm
    one_γr2 = 1 + _gauss_pre * sqmahab # (1 + γr²)

    # second partials w.r.t. position
    dxx = 2 * Ag * γ_f2 * (1 + 2 * _gauss_pre * dx^2 / fwhm²)
    dyy = 2 * Ag * γ_f2 * (1 + 2 * _gauss_pre * dy^2 / fwhm²)
    dxy = 4 * Ag * _gauss_pre^2 * dx * dy / fwhm²^2

    # position × fwhm
    dxf = 4 * Ag * _gauss_pre * dx / fwhm^3 * (2 + _gauss_pre * sqmahab)
    dyf = 4 * Ag * _gauss_pre * dy / fwhm^3 * (2 + _gauss_pre * sqmahab)

    # position × flux
    dxfl = -2 * _gauss_pre * g * dx / (fwhm² * norm)
    dyfl = -2 * _gauss_pre * g * dy / (fwhm² * norm)

    # fwhm × fwhm
    dff  = Ag / fwhm² * (4 * one_γr2^2 + 2 * (1 + 3 * _gauss_pre * sqmahab))

    # fwhm × flux
    dff_l = -2 * g * one_γr2 / (fwhm * norm)

    # assemble symmetric 5×5 matrix
    H = SA[
        dxx   dxy   dxf    dxfl   0
        dxy   dyy   dyf    dyfl   0
        dxf   dyf   dff    dff_l  0
        dxfl  dyfl  dff_l  0      0
        0     0     0      0      0
    ]
    return H
end

############################################

"""
    GaussianPSF(x, y, x_fwhm, y_fwhm, theta, flux, bkg) -> GaussianPSF{T}

General asymmetric Gaussian PSF with centroid `(x, y)`, FWHM along the x and y axes given by `x_fwhm` and `y_fwhm`, and rotated by `theta` degrees counter-clockwise from the x-axis. The `flux` is the integral of the PSF over all space, and `bkg` is a scalar background level added to the PSF.

```jldoctest
julia> using PSFModels: GaussianPSF

julia> GaussianPSF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=4.0, bkg=5.0) isa GaussianPSF{Float64}
true
```
"""
Base.@kwdef struct GaussianPSF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    x_fwhm::T
    y_fwhm::T
    theta::T
    flux::T
    bkg::T
    function GaussianPSF(x, y, x_fwhm, y_fwhm, theta, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(x_fwhm), typeof(y_fwhm), typeof(theta), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T # promote to Float if all inputs are integers
        return new{T}(T(x), T(y), T(x_fwhm), T(y_fwhm), T(theta), T(flux), T(bkg))
    end
end
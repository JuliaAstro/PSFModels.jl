
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

# import ForwardDiff: gradient
# using PSFModels: CircularGaussianPSF, fit_deriv, evaluate
# using ConstructionBase: getproperties
# t = CircularGaussianPSF(x=1.0, y=2.0, fwhm=3.0, flux=6.0, bkg=7.0)
# gradient(x->evaluate(CircularGaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(fit_deriv(t, -1, 1))
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

# import ForwardDiff: hessian
# using PSFModels: CircularGaussianPSF, fit_hessian, evaluate
# using ConstructionBase: getproperties
# t = CircularGaussianPSF(x=1.0, y=2.0, fwhm=3.0, flux=6.0, bkg=7.0)
# hessian(x->evaluate(CircularGaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(fit_hessian(t, -1, 1))
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

function evaluate(model::GaussianPSF{T}, px, py) where T
    d = T(π) / 180
    dx = px - model.x
    dy = py - model.y
    cs = cos(d * model.theta)
    sn = sin(d * model.theta)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    # u, v = rotate_point(px - model.x, py - model.y, model.theta)
    ax = model.x_fwhm
    ay = model.y_fwhm
    sqmahab = (u / ax)^2 + (v / ay)^2
    norm = -(T(π) * ax * ay / T(GAUSS_PRE))
    amp = model.flux / norm
    return muladd(amp, exp(T(GAUSS_PRE) * sqmahab), model.bkg)
end

# import ForwardDiff: gradient
# using PSFModels: GaussianPSF, fit_deriv, evaluate
# using ConstructionBase: getproperties
# t = GaussianPSF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=6.0, bkg=7.0)
# gradient(x->evaluate(GaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(fit_deriv(t, -1, 1))
function fit_deriv(model::GaussianPSF{T}, px, py) where T
    _gauss_pre = T(GAUSS_PRE)
    d = T(π) / 180
    dx = px - model.x
    dy = py - model.y
    ax = model.x_fwhm
    ay = model.y_fwhm
    ax² = ax^2
    ay² = ay^2
    cs = cos(d * model.theta)
    sn = sin(d * model.theta)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    sqmahab = u^2 / ax² + v^2 / ay²
    norm = -(T(π) * ax * ay / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)
    Ag = amp * g
    # ∂f/∂x, ∂f/∂y
    df_dx     = -2 * Ag * _gauss_pre * (cs * u / ax² - sn * v / ay²)
    df_dy     = -2 * Ag * _gauss_pre * (sn * u / ax² + cs * v / ay²)
    # ∂f/∂x_fwhm, ∂f/∂y_fwhm
    df_dax    = -Ag / ax * (1 + 2 * _gauss_pre * u^2 / ax²)
    df_day    = -Ag / ay * (1 + 2 * _gauss_pre * v^2 / ay²)
    # ∂f/∂theta (theta in degrees; d = π/180 factor)
    df_dtheta = 2 * d * Ag * _gauss_pre * u * v * (1 / ax² - 1 / ay²)
    # ∂f/∂flux
    df_dflux  = g / norm
    # ∂f/∂bkg
    df_dbkg   = one(T)
    return SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]
end

# import ForwardDiff: hessian
# using PSFModels: GaussianPSF, fit_hessian, evaluate
# using ConstructionBase: getproperties
# t = GaussianPSF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=6.0, bkg=7.0)
# hessian(x->evaluate(GaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(fit_hessian(t, -1, 1))
function fit_hessian(model::GaussianPSF{T}, px, py) where T
    _gauss_pre = T(GAUSS_PRE)
    d = T(π) / 180
    dx = px - model.x
    dy = py - model.y
    ax = model.x_fwhm
    ay = model.y_fwhm
    ax² = ax^2
    ay² = ay^2
    cs = cos(d * model.theta)
    sn = sin(d * model.theta)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    sqmahab = u^2 / ax² + v^2 / ay²
    norm = -(T(π) * ax * ay / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)
    Ag = amp * g
    γ = _gauss_pre
    D = 1 / ax² - 1 / ay²  # (1/ax² - 1/ay²), used in theta terms

    # First derivatives of sqmahab w.r.t. each parameter
    Qx     = -2 * (cs * u / ax² - sn * v / ay²)
    Qy     = -2 * (sn * u / ax² + cs * v / ay²)
    Qax    = -2 * u^2 / ax^3
    Qay    = -2 * v^2 / ay^3
    Qtheta = d * 2 * u * v * D

    # Second derivatives of sqmahab
    Rxx      = 2 * (cs^2 / ax² + sn^2 / ay²)
    Ryy      = 2 * (sn^2 / ax² + cs^2 / ay²)
    Rxy      = 2 * cs * sn * D
    Rxax     = 4 * cs * u / ax^3
    Rxay     = -4 * sn * v / ay^3
    Ryax     = 4 * sn * u / ax^3
    Ryay     = 4 * cs * v / ay^3
    Raxax    = 6 * u^2 / ax^4
    Rayay    = 6 * v^2 / ay^4
    Rxtheta  = d * 2 * (sn * u - cs * v) * D
    Rytheta  = -d * 2 * (cs * u + sn * v) * D
    Raxtheta = -d * 4 * u * v / ax^3
    Raytheta = d * 4 * u * v / ay^3
    Rtheta2  = d^2 * 2 * (v^2 - u^2) * D

    # Hessian entries (row/col order: x, y, x_fwhm, y_fwhm, theta, flux, bkg)
    # Second derivatives involving only sqmahab (∂A/∂x = ∂A/∂y = ∂A/∂theta = 0)
    dxx      = Ag * (γ^2 * Qx^2     + γ * Rxx)
    dyy      = Ag * (γ^2 * Qy^2     + γ * Ryy)
    dxy      = Ag * (γ^2 * Qx * Qy  + γ * Rxy)
    dxtheta  = Ag * (γ^2 * Qx * Qtheta  + γ * Rxtheta)
    dytheta  = Ag * (γ^2 * Qy * Qtheta  + γ * Rytheta)
    dtheta2  = Ag * (γ^2 * Qtheta^2 + γ * Rtheta2)

    # Cross terms: position × fwhm (∂A/∂ax = -A/ax, ∂A/∂ay = -A/ay)
    dxax     = Ag * γ * (γ * Qx * Qax + Rxax  - Qx / ax)
    dxay     = Ag * γ * (γ * Qx * Qay + Rxay  - Qx / ay)
    dyax     = Ag * γ * (γ * Qy * Qax + Ryax  - Qy / ax)
    dyay     = Ag * γ * (γ * Qy * Qay + Ryay  - Qy / ay)

    # Cross terms: fwhm × theta (∂A/∂theta = 0)
    daxtheta = Ag * γ * (γ * Qax * Qtheta + Raxtheta - Qtheta / ax)
    daytheta = Ag * γ * (γ * Qay * Qtheta + Raytheta - Qtheta / ay)

    # fwhm × fwhm (∂²A/∂ax² = 2A/ax², ∂²A/∂ay² = 2A/ay², ∂²A/∂ax∂ay = A/(ax*ay))
    daxax    = Ag * (2 / ax² + γ^2 * Qax^2 + γ * Raxax - 2 * γ * Qax / ax)
    dayay    = Ag * (2 / ay² + γ^2 * Qay^2 + γ * Rayay - 2 * γ * Qay / ay)
    daxay    = Ag * (1 / (ax * ay) - γ * Qay / ax - γ * Qax / ay + γ^2 * Qax * Qay)

    # Cross terms with flux (∂A/∂flux = 1/norm, ∂²A/∂flux∂ax = 1/(norm*ax), etc.)
    dxflux   = γ * g * Qx / norm
    dyflux   = γ * g * Qy / norm
    daxflux  = g / norm * (-1 / ax + γ * Qax)
    dayflux  = g / norm * (-1 / ay + γ * Qay)
    dthetaflux = γ * g * Qtheta / norm

    H = SA[
        dxx     dxy     dxax     dxay     dxtheta     dxflux     0
        dxy     dyy     dyax     dyay     dytheta     dyflux     0
        dxax    dyax    daxax    daxay    daxtheta    daxflux    0
        dxay    dyay    daxay    dayay    daytheta    dayflux    0
        dxtheta dytheta daxtheta daytheta dtheta2     dthetaflux 0
        dxflux  dyflux  daxflux  dayflux  dthetaflux  0          0
        0       0       0        0        0           0          0
    ]
    return H
end
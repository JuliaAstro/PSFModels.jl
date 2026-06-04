"""Gaussian exponent coefficient when parameterized by FWHM:
exp(-4*log(2) * r² / fwhm²)  = exp(GAUSS_PRE * r² / fwhm²) for effective radius r"""
const GAUSS_PRE = -4 * log(2)

"""Solution to ``J_1(π R_z) = 0`` where ``J_1`` is the first-order Bessel function of the first kind. This is used to compute the effective radius of an Airy disk in terms of its FWHM."""
const AIRY_RZ = 1.2196698912665045

@doc raw"""
    CircularGaussianPSF(x, y, fwhm, flux, bkg) → CircularGaussianPSF{T}

Circular, symmetric Gaussian PSF with centroid `(x, y)` and FWHM given by `fwhm`. 
The `flux` is the integral of the PSF over all space, and `bkg` is a scalar
background level added to the PSF. The model is evaluated by sampling the 2D
circularly-symmetric Gaussian function at the given position and adding the background:

```math
I(x, y) = \frac{F}{\pi\,\mathrm{FWHM}^2/(4\ln 2)}
\exp\!\left[-4\ln 2\,
\frac{(x-x_0)^2 + (y-y_0)^2}{\mathrm{FWHM}^2}\right]
+ B
```

where ``(x_0, y_0)`` is the centroid, ``F`` is the total flux, and ``B`` is the background level.

# Examples
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
peak(model::CircularGaussianPSF{T}) where {T} = model.flux / (π * model.fwhm^2 / -T(GAUSS_PRE)) + model.bkg
effective_area(model::CircularGaussianPSF{T}) where {T} = π * model.fwhm^2 / T(2 * log(2))

function evaluate(model::CircularGaussianPSF{T}, px, py) where {T}
    dx = px - model.x
    dy = py - model.y
    sqmahab = (dx^2 + dy^2) / model.fwhm^2
    norm = -(π * model.fwhm^2 / T(GAUSS_PRE))
    amp = model.flux / norm
    return muladd(amp, exp(T(GAUSS_PRE) * sqmahab), model.bkg)
end

# using ForwardDiff: gradient
# using PSFModels: CircularGaussianPSF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = CircularGaussianPSF(x=1.0, y=2.0, fwhm=3.0, flux=6.0, bkg=7.0)
# _, g = evaluate_fg(t, -1, 1)
# gradient(x->evaluate(CircularGaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::CircularGaussianPSF{T}, px, py) where {T}
    _gauss_pre = T(GAUSS_PRE)
    dx = px - model.x
    dy = py - model.y
    fwhm = model.fwhm
    fwhm² = fwhm^2
    sqmahab = (dx^2 + dy^2) / fwhm²
    norm = -(T(π) * fwhm² / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)
    Ag = amp * g
    f = muladd(amp, g, model.bkg)
    # Gradient ↓
    γ_f2 = _gauss_pre / fwhm²
    df_dx = -2 * Ag * γ_f2 * dx
    df_dy = -2 * Ag * γ_f2 * dy
    df_dfwhm = -2 * Ag * (1 + _gauss_pre * sqmahab) / fwhm
    df_dflux = g / norm
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dfwhm, df_dflux, df_dbkg]
    return f, G
end

# using ForwardDiff: hessian
# using PSFModels: CircularGaussianPSF, evaluate_fgh, evaluate
# using ConstructionBase: getproperties
# t = CircularGaussianPSF(x=1.0, y=2.0, fwhm=3.0, flux=6.0, bkg=7.0)
# _, _, h = evaluate_fgh(t, -1, 1)
# hessian(x->evaluate(CircularGaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(h)
function evaluate_fgh(model::CircularGaussianPSF{T}, px, py) where {T}
    _gauss_pre = T(GAUSS_PRE)
    dx = px - model.x
    dy = py - model.y
    fwhm = model.fwhm
    fwhm² = fwhm^2
    sqmahab = (dx^2 + dy^2) / fwhm²
    norm = -(T(π) * fwhm² / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)
    Ag = amp * g
    f = muladd(amp, g, model.bkg)

    # Shared sub-expressions
    γ_f2 = _gauss_pre / fwhm²  # γ/fwhm²
    one_γr2 = 1 + _gauss_pre * sqmahab  # (1 + γr²)

    # Gradient
    df_dx = -2 * Ag * γ_f2 * dx
    df_dy = -2 * Ag * γ_f2 * dy
    df_dfwhm = -2 * Ag * one_γr2 / fwhm
    df_dflux = g / norm
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dfwhm, df_dflux, df_dbkg]

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
    dff = Ag / fwhm² * (4 * one_γr2^2 + 2 * (1 + 3 * _gauss_pre * sqmahab))

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
    return f, G, H
end

############################################

@doc raw"""
    GaussianPSF(x, y, x_fwhm, y_fwhm, theta, flux, bkg) -> GaussianPSF{T}

General asymmetric Gaussian PSF with centroid `(x, y)`, FWHM along the x and y axes
given by `x_fwhm` and `y_fwhm`, and rotated by `theta` degrees counter-clockwise
from the x-axis. The `flux` is the integral of the PSF over all space, and `bkg`
is a scalar background level added to the PSF.

The model is evaluated by sampling the general 2D
Gaussian function at the given position and adding the background:

```math
I(x, y) = \frac{F}{\pi\,\mathrm{FWHM}_x\,\mathrm{FWHM}_y/(4\ln 2)}
\exp\!\left[-4\ln 2\,
\left(
\frac{u^2}{\mathrm{FWHM}_x^2}
+
\frac{v^2}{\mathrm{FWHM}_y^2}
\right)\right]
+ B,
```

where

```math
u = \cos\theta\,(x-x_0) + \sin\theta\,(y-y_0),
\qquad
v = -\sin\theta\,(x-x_0) + \cos\theta\,(y-y_0),
```

``F`` is the total flux and ``B`` is the background level.

# Examples
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
peak(model::GaussianPSF{T}) where {T} = model.flux / (π * model.x_fwhm * model.y_fwhm / -T(GAUSS_PRE)) + model.bkg
effective_area(model::GaussianPSF) = π * model.x_fwhm * model.y_fwhm / (2 * log(2))

function evaluate(model::GaussianPSF{T}, px, py) where {T}
    θ = deg2rad(model.theta)
    dx = px - model.x
    dy = py - model.y
    sn, cs = sincos(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    ax = model.x_fwhm
    ay = model.y_fwhm
    sqmahab = (u / ax)^2 + (v / ay)^2
    norm = -(T(π) * ax * ay / T(GAUSS_PRE))
    amp = model.flux / norm
    return muladd(amp, exp(T(GAUSS_PRE) * sqmahab), model.bkg)
end

# using ForwardDiff: gradient
# using PSFModels: GaussianPSF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = GaussianPSF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=6.0, bkg=7.0)
# _, g = evaluate_fg(t, -1, 1)
# gradient(x->evaluate(GaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::GaussianPSF{T}, px, py) where {T}
    _gauss_pre = T(GAUSS_PRE)
    d = deg2rad(one(T))  # factor to convert degrees to radians for theta derivatives
    dx = px - model.x
    dy = py - model.y
    ax = model.x_fwhm
    ay = model.y_fwhm
    ax² = ax^2
    ay² = ay^2
    θ = deg2rad(model.theta)
    sn, cs = sincos(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    sqmahab = u^2 / ax² + v^2 / ay²
    norm = -(T(π) * ax * ay / _gauss_pre)
    amp = model.flux / norm
    g = exp(_gauss_pre * sqmahab)
    Ag = amp * g
    f = muladd(amp, g, model.bkg)

    # Gradient ↓
    γ = _gauss_pre
    D = 1 / ax² - 1 / ay²
    Qx = -2 * (cs * u / ax² - sn * v / ay²)
    Qy = -2 * (sn * u / ax² + cs * v / ay²)
    Qax = -2 * u^2 / ax^3
    Qay = -2 * v^2 / ay^3
    Qtheta = d * 2 * u * v * D
    df_dx = Ag * γ * Qx
    df_dy = Ag * γ * Qy
    df_dax = Ag * (-1 / ax + γ * Qax)
    df_day = Ag * (-1 / ay + γ * Qay)
    df_dtheta = Ag * γ * Qtheta
    df_dflux = g / norm
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]
    return f, G
end

# using ForwardDiff: hessian
# using PSFModels: GaussianPSF, evaluate_fgh, evaluate
# using ConstructionBase: getproperties
# t = GaussianPSF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=6.0, bkg=7.0)
# _, _, h = evaluate_fgh(t, -1, 1)
# hessian(x->evaluate(GaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(h)
function evaluate_fgh(model::GaussianPSF{T}, px, py) where {T}
    γ = T(GAUSS_PRE)
    d = deg2rad(one(T))  # factor to convert degrees to radians for theta derivatives
    dx = px - model.x
    dy = py - model.y
    ax = model.x_fwhm
    ay = model.y_fwhm
    ax² = ax^2
    ay² = ay^2
    θ = deg2rad(model.theta)
    cs = cos(θ)
    sn = sin(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    sqmahab = u^2 / ax² + v^2 / ay²
    norm = -(T(π) * ax * ay / γ)
    amp = model.flux / norm
    g = exp(γ * sqmahab)
    Ag = amp * g
    f = muladd(amp, g, model.bkg)

    # First derivatives of sqmahab w.r.t. each parameter (used for both gradient and Hessian)
    D = 1 / ax² - 1 / ay²  # (1/ax² - 1/ay²), used in theta terms
    Qx = -2 * (cs * u / ax² - sn * v / ay²)
    Qy = -2 * (sn * u / ax² + cs * v / ay²)
    Qax = -2 * u^2 / ax^3
    Qay = -2 * v^2 / ay^3
    Qtheta = d * 2 * u * v * D

    # Gradient
    df_dx = Ag * γ * Qx
    df_dy = Ag * γ * Qy
    df_dax = Ag * (-1 / ax + γ * Qax)
    df_day = Ag * (-1 / ay + γ * Qay)
    df_dtheta = Ag * γ * Qtheta
    df_dflux = g / norm
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]

    # Second derivatives of sqmahab
    Rxx = 2 * (cs^2 / ax² + sn^2 / ay²)
    Ryy = 2 * (sn^2 / ax² + cs^2 / ay²)
    Rxy = 2 * cs * sn * D
    Rxax = 4 * cs * u / ax^3
    Rxay = -4 * sn * v / ay^3
    Ryax = 4 * sn * u / ax^3
    Ryay = 4 * cs * v / ay^3
    Raxax = 6 * u^2 / ax^4
    Rayay = 6 * v^2 / ay^4
    Rxtheta = d * 2 * (sn * u - cs * v) * D
    Rytheta = -d * 2 * (cs * u + sn * v) * D
    Raxtheta = -d * 4 * u * v / ax^3
    Raytheta = d * 4 * u * v / ay^3
    Rtheta2 = d^2 * 2 * (v^2 - u^2) * D

    # Hessian entries (row/col order: x, y, x_fwhm, y_fwhm, theta, flux, bkg)
    # Second derivatives involving only sqmahab (∂A/∂x = ∂A/∂y = ∂A/∂theta = 0)
    dxx = Ag * (γ^2 * Qx^2 + γ * Rxx)
    dyy = Ag * (γ^2 * Qy^2 + γ * Ryy)
    dxy = Ag * (γ^2 * Qx * Qy + γ * Rxy)
    dxtheta = Ag * (γ^2 * Qx * Qtheta + γ * Rxtheta)
    dytheta = Ag * (γ^2 * Qy * Qtheta + γ * Rytheta)
    dtheta2 = Ag * (γ^2 * Qtheta^2 + γ * Rtheta2)

    # Cross terms: position × fwhm (∂A/∂ax = -A/ax, ∂A/∂ay = -A/ay)
    dxax = Ag * γ * (γ * Qx * Qax + Rxax - Qx / ax)
    dxay = Ag * γ * (γ * Qx * Qay + Rxay - Qx / ay)
    dyax = Ag * γ * (γ * Qy * Qax + Ryax - Qy / ax)
    dyay = Ag * γ * (γ * Qy * Qay + Ryay - Qy / ay)

    # Cross terms: fwhm × theta (∂A/∂theta = 0)
    daxtheta = Ag * γ * (γ * Qax * Qtheta + Raxtheta - Qtheta / ax)
    daytheta = Ag * γ * (γ * Qay * Qtheta + Raytheta - Qtheta / ay)

    # fwhm × fwhm (∂²A/∂ax² = 2A/ax², ∂²A/∂ay² = 2A/ay², ∂²A/∂ax∂ay = A/(ax*ay))
    daxax = Ag * (2 / ax² + γ^2 * Qax^2 + γ * Raxax - 2 * γ * Qax / ax)
    dayay = Ag * (2 / ay² + γ^2 * Qay^2 + γ * Rayay - 2 * γ * Qay / ay)
    daxay = Ag * (1 / (ax * ay) - γ * Qay / ax - γ * Qax / ay + γ^2 * Qax * Qay)

    # Cross terms with flux (∂A/∂flux = 1/norm, ∂²A/∂flux∂ax = 1/(norm*ax), etc.)
    dxflux = γ * g * Qx / norm
    dyflux = γ * g * Qy / norm
    daxflux = g / norm * (-1 / ax + γ * Qax)
    dayflux = g / norm * (-1 / ay + γ * Qay)
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
    return f, G, H
end

############################################

@doc raw"""
    CircularGaussianPRF(x, y, fwhm, flux, bkg) -> CircularGaussianPRF{T}

Circular, symmetric Gaussian pixel response function (PRF) with centroid `(x, y)` and FWHM
given by `fwhm`. The PRF is the underlying Gaussian PSF integrated analytically over each
pixel. The `flux` is the total flux (sum of PRF values over all pixels equals `flux`), and
`bkg` is a scalar background level added to the PRF.

The PRF value at pixel center `(px, py)` is the integral of the Gaussian PSF over the pixel
area `[px-0.5, px+0.5] × [py-0.5, py+0.5]`. The model is evaluated as

```math
f(px, py) = \frac{\mathrm{flux}}{4}
    \left[\mathrm{erf}\!\left(\frac{2\sqrt{\ln 2}\,(px + 0.5 - x)}{\mathrm{fwhm}}\right)
         -\mathrm{erf}\!\left(\frac{2\sqrt{\ln 2}\,(px - 0.5 - x)}{\mathrm{fwhm}}\right)\right]
    \left[\mathrm{erf}\!\left(\frac{2\sqrt{\ln 2}\,(py + 0.5 - y)}{\mathrm{fwhm}}\right)
         -\mathrm{erf}\!\left(\frac{2\sqrt{\ln 2}\,(py - 0.5 - y)}{\mathrm{fwhm}}\right)\right]
    + \mathrm{bkg}
```

# Examples
```jldoctest
julia> using PSFModels: CircularGaussianPRF

julia> CircularGaussianPRF(1.0, 2.0, 3.0, 4.0, 5.0) isa CircularGaussianPRF{Float64}
true
```
"""
Base.@kwdef struct CircularGaussianPRF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    fwhm::T
    flux::T
    bkg::T
    function CircularGaussianPRF(x, y, fwhm, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(fwhm), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T
        return new{T}(T(x), T(y), T(fwhm), T(flux), T(bkg))
    end
end

# Peak occurs when centroid is exactly at a pixel centre: flux * erf(√ln2 / fwhm)² + bkg
peak(model::CircularGaussianPRF{T}) where {T} =
    model.flux * erf(sqrt(T(log(2))) / model.fwhm)^2 + model.bkg
effective_area(model::CircularGaussianPRF{T}) where {T} = π * model.fwhm^2 / T(2 * log(2))

function evaluate(model::CircularGaussianPRF{T}, px, py) where {T}
    α = 2 * sqrt(T(log(2))) / model.fwhm
    dx = px - model.x
    dy = py - model.y
    Ex = erf(α * (dx + T(0.5))) - erf(α * (dx - T(0.5)))
    Ey = erf(α * (dy + T(0.5))) - erf(α * (dy - T(0.5)))
    return muladd(model.flux / 4, Ex * Ey, model.bkg)
end

# using ForwardDiff: gradient
# using PSFModels: CircularGaussianPRF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = CircularGaussianPRF(x=0.0, y=0.0, fwhm=10.0, flux=1.0, bkg=10.0)
# _, g = evaluate_fg(t, 1, 2)
# gradient(x->evaluate(CircularGaussianPRF(x...), 1, 2), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::CircularGaussianPRF{T}, px, py) where {T}
    α = 2 * sqrt(T(log(2))) / model.fwhm
    dx = px - model.x
    dy = py - model.y
    u_p = α * (dx + T(0.5))
    u_m = α * (dx - T(0.5))
    v_p = α * (dy + T(0.5))
    v_m = α * (dy - T(0.5))
    Ex = erf(u_p) - erf(u_m)
    Ey = erf(v_p) - erf(v_m)
    f = muladd(model.flux / 4, Ex * Ey, model.bkg)

    # Derivative factors: d/du erf(u) = (2/√π) exp(-u²)
    _2_sqrtpi = 2 / sqrt(T(π))
    Gxp = _2_sqrtpi * exp(-u_p^2)
    Gxm = _2_sqrtpi * exp(-u_m^2)
    Gyp = _2_sqrtpi * exp(-v_p^2)
    Gym = _2_sqrtpi * exp(-v_m^2)

    df_dx = model.flux / 4 * Ey * α * (Gxm - Gxp)
    df_dy = model.flux / 4 * Ex * α * (Gym - Gyp)
    df_dfwhm = model.flux / 4 / model.fwhm * ((Gxm * u_m - Gxp * u_p) * Ey + Ex * (Gym * v_m - Gyp * v_p))
    df_dflux = Ex * Ey / 4
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dfwhm, df_dflux, df_dbkg]
    return f, G
end

############################################

@doc raw"""
    GaussianPRF(x, y, x_fwhm, y_fwhm, theta, flux, bkg) -> GaussianPRF{T}

Asymmetric Gaussian pixel response function (PRF) with centroid `(x, y)`, FWHM along the
x and y axes given by `x_fwhm` and `y_fwhm`, and rotated by `theta` degrees
counter-clockwise from the x-axis. The PRF is the underlying Gaussian PSF
integrated analytically over each pixel. The `flux` is the total flux (sum of PRF values
over all pixels equals `flux`), and `bkg` is a scalar background level added to the PRF.

The PRF value at pixel center `(px, py)` is the integral of the Gaussian PSF over the pixel
area `[px-0.5, px+0.5] × [py-0.5, py+0.5]`, evaluated after rotating the coordinates
by `-theta` to align with the principal axes. The model is evaluated as 

```math
f(px, py) = \\frac{\\mathrm{flux}}{4}
    \\left[\\mathrm{erf}\\!\\left(\\frac{2\\sqrt{\\ln 2}\\,(u + 0.5)}{x\\_\\mathrm{fwhm}}\\right)
         -\\mathrm{erf}\\!\\left(\\frac{2\\sqrt{\\ln 2}\\,(u - 0.5)}{x\\_\\mathrm{fwhm}}\\right)\\right]
    \\left[\\mathrm{erf}\\!\\left(\\frac{2\\sqrt{\\ln 2}\\,(v + 0.5)}{y\\_\\mathrm{fwhm}}\\right)
         -\\mathrm{erf}\\!\\left(\\frac{2\\sqrt{\\ln 2}\\,(v - 0.5)}{y\\_\\mathrm{fwhm}}\\right)\\right]
    + \\mathrm{bkg}
```

where 

```math
u = \cos\theta\,(x-x_0) + \sin\theta\,(y-y_0),
\qquad
v = -\sin\theta\,(x-x_0) + \cos\theta\,(y-y_0),
```

```jldoctest
julia> using PSFModels: GaussianPRF

julia> GaussianPRF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=4.0, bkg=5.0) isa GaussianPRF{Float64}
true
```
"""
Base.@kwdef struct GaussianPRF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    x_fwhm::T
    y_fwhm::T
    theta::T
    flux::T
    bkg::T
    function GaussianPRF(x, y, x_fwhm, y_fwhm, theta, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(x_fwhm), typeof(y_fwhm), typeof(theta), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T
        return new{T}(T(x), T(y), T(x_fwhm), T(y_fwhm), T(theta), T(flux), T(bkg))
    end
end

# Peak occurs when centroid is at a pixel centre: flux * erf(√ln2/x_fwhm) * erf(√ln2/y_fwhm) + bkg
peak(model::GaussianPRF{T}) where {T} =
    model.flux * erf(sqrt(T(log(2))) / model.x_fwhm) * erf(sqrt(T(log(2))) / model.y_fwhm) + model.bkg
effective_area(model::GaussianPRF{T}) where {T} = π * model.x_fwhm * model.y_fwhm / T(2 * log(2))

function evaluate(model::GaussianPRF{T}, px, py) where {T}
    c = sqrt(-T(GAUSS_PRE)) # 2 * sqrt(T(log(2)))
    αx = c / model.x_fwhm
    αy = c / model.y_fwhm
    dx = px - model.x
    dy = py - model.y
    # rotate coordinates into axis-aligned frame
    θ = deg2rad(model.theta)
    sn, cs = sincos(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    Ex = erf(αx * (u + T(0.5))) - erf(αx * (u - T(0.5)))
    Ey = erf(αy * (v + T(0.5))) - erf(αy * (v - T(0.5)))
    return muladd(model.flux / 4, Ex * Ey, model.bkg)
end

# using ForwardDiff: gradient
# using PSFModels: GaussianPRF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = GaussianPRF(x=0.0, y=0.0, x_fwhm=10.0, y_fwhm=6.0, theta=45.0, flux=1.0, bkg=10.0)
# _, g = evaluate_fg(t, 1, 2)
# gradient(x->evaluate(GaussianPRF(x...), 1, 2), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::GaussianPRF{T}, px, py) where {T}
    c = sqrt(-T(GAUSS_PRE)) # 2 * sqrt(T(log(2)))
    αx = c / model.x_fwhm
    αy = c / model.y_fwhm
    dx = px - model.x
    dy = py - model.y
    # rotate coordinates into axis-aligned frame
    θ = deg2rad(model.theta)
    sn, cs = sincos(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    u_p = αx * (u + T(0.5))
    u_m = αx * (u - T(0.5))
    v_p = αy * (v + T(0.5))
    v_m = αy * (v - T(0.5))
    Ex = erf(u_p) - erf(u_m)
    Ey = erf(v_p) - erf(v_m)
    f = muladd(model.flux / 4, Ex * Ey, model.bkg)

    # Derivative factors: d/du erf(u) = (2/√π) exp(-u²)
    _2_sqrtpi = T(2 / sqrt(π))
    Gxp = _2_sqrtpi * exp(-u_p^2)
    Gxm = _2_sqrtpi * exp(-u_m^2)
    Gyp = _2_sqrtpi * exp(-v_p^2)
    Gym = _2_sqrtpi * exp(-v_m^2)

    fl4 = model.flux / 4
    # partials of the difference-of-erfs w.r.t. the rotated coordinates
    dEx_du = αx * (Gxm - Gxp)
    dEy_dv = αy * (Gym - Gyp)

    # ∂f/∂x: ∂u/∂x = -cs, ∂v/∂x = sn
    #   ∂f/∂x = fl4*((-dEx_du)*(-cs)*Ey + Ex*(-dEy_dv)*sn) = fl4*(dEx_du*cs*Ey - dEy_dv*sn*Ex)
    df_dx = fl4 * (cs * dEx_du * Ey - sn * dEy_dv * Ex)
    # ∂f/∂y: ∂u/∂y = -sn, ∂v/∂y = -cs
    #   ∂f/∂y = fl4*((-dEx_du)*(-sn)*Ey + Ex*(-dEy_dv)*(-cs)) = fl4*(dEx_du*sn*Ey + dEy_dv*cs*Ex)
    df_dy = fl4 * (sn * dEx_du * Ey + cs * dEy_dv * Ex)
    # ∂f/∂x_fwhm: same form but with u replacing dx
    df_dax = fl4 / model.x_fwhm * (Gxm * u_m - Gxp * u_p) * Ey
    df_day = fl4 / model.y_fwhm * Ex * (Gym * v_m - Gyp * v_p)
    # ∂f/∂θ: ∂u/∂θ = d*v, ∂v/∂θ = -d*u, where d = deg2rad(1) = π / 180
    #   ∂f/∂θ = fl4*((-dEx_du)*d*v*Ey + Ex*(-dEy_dv)*(-d*u)) = fl4*d*(dEy_dv*u*Ex - dEx_du*v*Ey)
    df_dtheta = fl4 * deg2rad(one(T)) * (dEy_dv * u * Ex - dEx_du * v * Ey)
    df_dflux = Ex * Ey / 4
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]
    return f, G
end

@doc raw"""
    AiryPSF(x, y, radius, flux, bkg) -> AiryPSF{T}

Airy disk point spread function (PSF) with centroid `(x, y)`, first-dark-ring `radius`,
total `flux`, and scalar background `bkg`.

The `radius` parameter is the radial distance from the centroid to the first dark ring
(also called the first zero) of the Airy pattern. This is the limiting angular
resolution of an optical system with ``R = R_z \lambda/D``, where `D` is the
aperture diameter, ``\lambda`` is the wavelength of the light, and `R_z ≈ 1.2197` is the
first positive zero of the first-order Bessel function `J_1`. 

The model is evaluated by sampling the Airy disk PSF
at the given position and adding the background:

```math
f(px, py) = \frac{\pi\,\mathrm{flux}}{4\,a^2}
            \left(\frac{2\,J_1(u)}{u}\right)^2 + \mathrm{bkg}
```

where

```math
r = \sqrt{(px - x)^2 + (py - y)^2}, \qquad
u = \frac{\pi\,r}{a}, \qquad
a = \frac{\mathrm{radius}}{r_z}, \qquad
r_z = \frac{j_{1,1}}{\pi} \approx 1.2197
```

with $j_{1,1}$ the first positive zero of the first-order Bessel function $J_1$.
The singularity at $u = 0$ is handled analytically; the peak value at the centroid
is $\pi\,\mathrm{flux}/(4a^2) + \mathrm{bkg}$.
The FWHM is approximately $0.8437 \times \mathrm{radius}$ along each axis.

```jldoctest
julia> using PSFModels: AiryPSF

julia> AiryPSF(x=1.0, y=2.0, radius=3.0, flux=4.0, bkg=5.0) isa AiryPSF{Float64}
true
```
"""
Base.@kwdef struct AiryPSF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    radius::T
    flux::T
    bkg::T

    function AiryPSF(x, y, radius, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(radius), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T
        return new{T}(T(x), T(y), T(radius), T(flux), T(bkg))
    end
end
amplitude(model::AiryPSF{T}) where {T} = model.flux / (model.radius / T(AIRY_RZ))^2 * T(π) / 4
peak(model::AiryPSF) = amplitude(model) + model.bkg
# 0.919... is 16 ∫ besselj1(u)^4 / u^3 du from 0 to ∞, numerically integrated
effective_area(model::AiryPSF{T}) where {T} = 8 * (model.radius / AIRY_RZ)^2 / π / T(0.9192407077670396)
function fwhm(model::AiryPSF{T}) where {T}
    _f = T(0.8436659602162364) * model.radius
    return (_f, _f)
end

function evaluate(model::AiryPSF{T}, px, py) where {T}
    # r = hypot(px - model.x, py - model.y) # hypot is slow...
    r = sqrt((px - model.x)^2 + (py - model.y)^2)
    a = model.radius / T(AIRY_RZ)
    u = π * r / a
    # Handle the u=0 case separately to avoid NaNs from besselj1(0)/0
    A2 = if abs(u) < eps(T)
        one(T)
    else
        J1 = besselj1(u)
        (2 * J1 / u)^2
    end
    norm = a^2 / π * 4
    amp = model.flux / norm
    return muladd(amp, A2, model.bkg)
end
# using ForwardDiff: gradient
# using PSFModels: AiryPSF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = AiryPSF(x=0.0, y=0.0, radius=10.0, flux=1.0, bkg=10.0)
# _, g = evaluate_fg(t, 1, 2)
# gradient(x->evaluate(AiryPSF(x...), 1, 2), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::AiryPSF{T}, px, py) where {T}
    dx = px - model.x
    dy = py - model.y
    # r = hypot(dx, dy) # hypot is slow...
    r = sqrt(dx^2 + dy^2)
    a = model.radius / T(AIRY_RZ)
    u = π * r / a
    if abs(u) < eps(T)
        A2 = one(T)
        dA2_du = zero(T)
    else
        J0 = besselj0(u)
        J1 = besselj1(u)
        J2 = besselj(2, u)
        A = 2J1 / u
        Ap = (u * (J0 - J2) - 2J1) / (u^2)
        A2 = A^2
        dA2_du = 2A * Ap
    end
    norm = a^2 / π * 4
    amp = model.flux / norm
    f = muladd(amp, A2, model.bkg)

    # Gradients ↓
    df_dflux = A2 / norm
    df_dbkg = one(T)
    # avoid division by zero at center
    if r == 0
        return f, SA[zero(T), zero(T), zero(T), df_dflux, df_dbkg]
    end
    du_dr = π / a
    df_dr = amp * dA2_du * du_dr
    df_dx = -df_dr * dx / r
    df_dy = -df_dr * dy / r
    da_dR = inv(T(AIRY_RZ))
    df_da = amp / a * (-u * dA2_du - 2 * A2)
    df_dradius = df_da * da_dR
    return f, SA[df_dx, df_dy, df_dradius, df_dflux, df_dbkg]
end

@doc raw"""
    CircularMoffatPSF(x, y, α, β, flux, bkg) -> CircularMoffatPSF{T}

Circular Moffat PSF with centroid `(x, y)`, scale length `α`, wing parameter
`β`, total `flux`, and background `bkg`.

The model is evaluated by sampling the Moffat PSF at the given position and
adding the background:

```math
I(x, y) =
F\,\frac{\beta - 1}{\pi\,\alpha^2}
\left(
1 + \frac{(x-x_0)^2 + (y-y_0)^2}{\alpha^2}
\right)^{-\beta}
+ B.
```

# Examples
```jldoctest
julia> using PSFModels: CircularMoffatPSF

julia> CircularMoffatPSF(x=1.0, y=2.0, α=3.0, β=2.5, flux=4.0, bkg=5.0) isa CircularMoffatPSF{Float64}
true
```
"""
Base.@kwdef struct CircularMoffatPSF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    α::T
    β::T
    flux::T
    bkg::T

    function CircularMoffatPSF(x, y, α, β, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(α), typeof(β), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T
        return new{T}(T(x), T(y), T(α), T(β), T(flux), T(bkg))
    end
end
amplitude(model::CircularMoffatPSF{T}) where {T} = model.flux * (model.β - 1) / (π * model.α^2)
peak(model::CircularMoffatPSF) = amplitude(model) + model.bkg
effective_area(model::CircularMoffatPSF{T}) where {T} = T(π) * model.α^2 * (2 * model.β - 1) / (model.β - 1)^2
function fwhm(model::CircularMoffatPSF{T}) where {T}
    _f = 2 * model.α * sqrt(exp2(1 / model.β) - 1)
    return (_f, _f)
end
function evaluate(model::CircularMoffatPSF{T}, px, py) where {T}
    r2 = (px - model.x)^2 + (py - model.y)^2
    amp = model.flux * (model.β - 1) / (π * model.α^2)
    return muladd(amp, (1 + r2 / model.α^2)^(-model.β), model.bkg)
end
# using ForwardDiff: gradient
# using PSFModels: CircularMoffatPSF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = CircularMoffatPSF(x=0.0, y=0.0, α=5.0, β=3.0, flux=50.0, bkg=10.0)
# _, g = evaluate_fg(t, 1, 2)
# gradient(x->evaluate(CircularMoffatPSF(x...), 1, 2), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::CircularMoffatPSF{T}, px, py) where {T}
    α, β = model.α, model.β
    α² = α^2
    dx = px - model.x
    dy = py - model.y
    r2 = dx^2 + dy^2
    u = 1 + r2 / α²
    profile = u^(-β)
    norm = T(π) * α² / (β - 1)
    amp = model.flux / norm
    Ag = amp * profile
    f = muladd(amp, profile, model.bkg)

    df_dx = 2 * Ag * β * dx / (α² * u)
    df_dy = 2 * Ag * β * dy / (α² * u)
    df_dα = Ag * (-2 / α + 2 * β * r2 / (α^3 * u))
    df_dβ = Ag * (1 / (β - 1) - log(u))
    df_dflux = profile / norm
    df_dbkg = one(T)
    return f, SA[df_dx, df_dy, df_dα, df_dβ, df_dflux, df_dbkg]
end

@doc raw"""
    MoffatPSF(x, y, x_α, y_α, theta, β, flux, bkg) -> MoffatPSF{T}

General asymmetric Moffat PSF with centroid `(x, y)`, scale lengths `x_α`
and `y_α` along the rotated axes, rotation angle `theta` in degrees
counter-clockwise from the x-axis, wing parameter `β`, total `flux`, and
background `bkg`.

The model is evaluated by sampling the Moffat PSF
at the given position and adding the background:

```math
I(x, y) =
F\,\frac{\beta - 1}{\pi\,\alpha_x\,\alpha_y}
\left(
1 + \frac{u^2}{\alpha_x^2} + \frac{v^2}{\alpha_y^2}
\right)^{-\beta}
+ B,
```

where

```math
u = \cos\theta\,(x-x_0) + \sin\theta\,(y-y_0),
\qquad
v = -\sin\theta\,(x-x_0) + \cos\theta\,(y-y_0).
```

# Examples
```jldoctest
julia> using PSFModels: MoffatPSF

julia> MoffatPSF(x=1.0, y=2.0, x_α=3.0, y_α=4.0, theta=35.0, β=2.5, flux=4.0, bkg=5.0) isa MoffatPSF{Float64}
true
```
"""
Base.@kwdef struct MoffatPSF{T} <: AbstractPSFModel{T}
    x::T
    y::T
    x_α::T
    y_α::T
    theta::T
    β::T
    flux::T
    bkg::T

    function MoffatPSF(x, y, x_α, y_α, theta, β, flux, bkg)
        T = promote_type(typeof(x), typeof(y), typeof(x_α), typeof(y_α), typeof(theta), typeof(β), typeof(flux), typeof(bkg))
        T = T <: Integer ? Float64 : T
        return new{T}(T(x), T(y), T(x_α), T(y_α), T(theta), T(β), T(flux), T(bkg))
    end
end
amplitude(model::MoffatPSF) = model.flux * (model.β - 1) / (π * model.x_α * model.y_α)
peak(model::MoffatPSF) = amplitude(model) + model.bkg
effective_area(model::MoffatPSF) = π * model.x_α * model.y_α * (2 * model.β - 1) / (model.β - 1)^2
function fwhm(model::MoffatPSF)
    f = 2 * sqrt(exp2(1 / model.β) - 1)
    return (model.x_α * f, model.y_α * f)
end
function evaluate(model::MoffatPSF{T}, px, py) where {T}
    θ = deg2rad(model.theta)
    dx = px - model.x
    dy = py - model.y
    sn, cs = sincos(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    ax = model.x_α
    ay = model.y_α
    q = u^2 / ax^2 + v^2 / ay^2
    amp = model.flux * (model.β - 1) / (T(π) * ax * ay)
    return muladd(amp, (1 + q)^(-model.β), model.bkg)
end
# using ForwardDiff: gradient
# using PSFModels: MoffatPSF, evaluate_fg, evaluate
# using ConstructionBase: getproperties
# t = MoffatPSF(x=1.0, y=2.0, x_α=3.0, y_α=4.0, theta=35.0, β=2.5, flux=6.0, bkg=7.0)
# _, g = evaluate_fg(t, -1, 1)
# gradient(x->evaluate(MoffatPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(g)
function evaluate_fg(model::MoffatPSF{T}, px, py) where {T}
    d = deg2rad(one(T))
    dx = px - model.x
    dy = py - model.y
    ax, ay, β = model.x_α, model.y_α, model.β
    ax² = ax^2
    ay² = ay^2
    θ = deg2rad(model.theta)
    sn, cs = sincos(θ)
    u = cs * dx + sn * dy
    v = -sn * dx + cs * dy
    q = u^2 / ax² + v^2 / ay²
    h = 1 + q
    profile = h^(-β)
    norm = π * ax * ay / (β - 1)
    amp = model.flux / norm
    Ag = amp * profile
    f = muladd(amp, profile, model.bkg)

    D = 1 / ax² - 1 / ay²
    Qx = -2 * (cs * u / ax² - sn * v / ay²)
    Qy = -2 * (sn * u / ax² + cs * v / ay²)
    Qax = -2 * u^2 / ax^3
    Qay = -2 * v^2 / ay^3
    Qtheta = d * 2 * u * v * D
    scale = -β / h
    df_dx = Ag * scale * Qx
    df_dy = Ag * scale * Qy
    df_dax = Ag * (-1 / ax + scale * Qax)
    df_day = Ag * (-1 / ay + scale * Qay)
    df_dtheta = Ag * scale * Qtheta
    df_dβ = Ag * (1 / (β - 1) - log(h))
    df_dflux = profile / norm
    df_dbkg = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dβ, df_dflux, df_dbkg]
    return f, G
end

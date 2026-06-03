# Gaussian exponent coefficient when parameterized by FWHM:
# exp(-4*log(2) * r² / fwhm²) for effective radius r
const GAUSS_PRE = -4 * log(2)

# factor for scaling radius in terms of the fwhm
const AIRY_PRE = π / 0.973

@doc raw"""
    CircularGaussianPSF(x, y, fwhm, flux, bkg) → CircularGaussianPSF{T}

Circular, symmetric Gaussian PSF with centroid `(x, y)` and FWHM given by `fwhm`. 
The `flux` is the integral of the PSF over all space, and `bkg` is a scalar background level added to the PSF.
The model is evaluated as 

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
peak(model::CircularGaussianPSF{T}) where T = model.flux / (π * model.fwhm^2 / -T(GAUSS_PRE)) + model.bkg
effective_area(model::CircularGaussianPSF{T}) where T = π * model.fwhm^2 / T(2 * log(2))

function evaluate(model::CircularGaussianPSF{T}, px, py) where T
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
function evaluate_fg(model::CircularGaussianPSF{T}, px, py) where T
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
    df_dx    = -2 * Ag * γ_f2 * dx
    df_dy    = -2 * Ag * γ_f2 * dy
    df_dfwhm = -2 * Ag * (1 + _gauss_pre * sqmahab) / fwhm
    df_dflux = g / norm
    df_dbkg  = one(T)
    G = SA[df_dx, df_dy, df_dfwhm, df_dflux, df_dbkg]
    return f, G
end

# using ForwardDiff: hessian
# using PSFModels: CircularGaussianPSF, evaluate_fgh, evaluate
# using ConstructionBase: getproperties
# t = CircularGaussianPSF(x=1.0, y=2.0, fwhm=3.0, flux=6.0, bkg=7.0)
# _, _, h = evaluate_fgh(t, -1, 1)
# hessian(x->evaluate(CircularGaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(h)
function evaluate_fgh(model::CircularGaussianPSF{T}, px, py) where T
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
    df_dx    = -2 * Ag * γ_f2 * dx
    df_dy    = -2 * Ag * γ_f2 * dy
    df_dfwhm = -2 * Ag * one_γr2 / fwhm
    df_dflux = g / norm
    df_dbkg  = one(T)
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
    return f, G, H
end

############################################

@doc raw"""
    GaussianPSF(x, y, x_fwhm, y_fwhm, theta, flux, bkg) -> GaussianPSF{T}

General asymmetric Gaussian PSF with centroid `(x, y)`, FWHM along the x and y axes
given by `x_fwhm` and `y_fwhm`, and rotated by `theta` degrees counter-clockwise
from the x-axis. The `flux` is the integral of the PSF over all space, and `bkg`
is a scalar background level added to the PSF.

The model is evaluated as

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
peak(model::GaussianPSF{T}) where T = model.flux / (π * model.x_fwhm * model.y_fwhm / -T(GAUSS_PRE)) + model.bkg
effective_area(model::GaussianPSF) = π * model.x_fwhm * model.y_fwhm / (2 * log(2))

function evaluate(model::GaussianPSF{T}, px, py) where T
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
function evaluate_fg(model::GaussianPSF{T}, px, py) where T
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
    Qx     = -2 * (cs * u / ax² - sn * v / ay²)
    Qy     = -2 * (sn * u / ax² + cs * v / ay²)
    Qax    = -2 * u^2 / ax^3
    Qay    = -2 * v^2 / ay^3
    Qtheta = d * 2 * u * v * D
    df_dx     = Ag * γ * Qx
    df_dy     = Ag * γ * Qy
    df_dax    = Ag * (-1 / ax + γ * Qax)
    df_day    = Ag * (-1 / ay + γ * Qay)
    df_dtheta = Ag * γ * Qtheta
    df_dflux  = g / norm
    df_dbkg   = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]
    return f, G
end

# using ForwardDiff: hessian
# using PSFModels: GaussianPSF, evaluate_fgh, evaluate
# using ConstructionBase: getproperties
# t = GaussianPSF(x=1.0, y=2.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=6.0, bkg=7.0)
# _, _, h = evaluate_fgh(t, -1, 1)
# hessian(x->evaluate(GaussianPSF(x...), -1, 1), collect(getproperties(t))) ≈ collect(h)
function evaluate_fgh(model::GaussianPSF{T}, px, py) where T
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
    Qx     = -2 * (cs * u / ax² - sn * v / ay²)
    Qy     = -2 * (sn * u / ax² + cs * v / ay²)
    Qax    = -2 * u^2 / ax^3
    Qay    = -2 * v^2 / ay^3
    Qtheta = d * 2 * u * v * D

    # Gradient
    df_dx     = Ag * γ * Qx
    df_dy     = Ag * γ * Qy
    df_dax    = Ag * (-1 / ax + γ * Qax)
    df_day    = Ag * (-1 / ay + γ * Qay)
    df_dtheta = Ag * γ * Qtheta
    df_dflux  = g / norm
    df_dbkg   = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]

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
    dxflux     = γ * g * Qx / norm
    dyflux     = γ * g * Qy / norm
    daxflux    = g / norm * (-1 / ax + γ * Qax)
    dayflux    = g / norm * (-1 / ay + γ * Qay)
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
area `[px-0.5, px+0.5] × [py-0.5, py+0.5]`.

The model is evaluated as 

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
peak(model::CircularGaussianPRF{T}) where T =
    model.flux * erf(sqrt(T(log(2))) / model.fwhm)^2 + model.bkg
effective_area(model::CircularGaussianPRF{T}) where T = π * model.fwhm^2 / T(2 * log(2))

function evaluate(model::CircularGaussianPRF{T}, px, py) where T
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
function evaluate_fg(model::CircularGaussianPRF{T}, px, py) where T
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
    df_dbkg  = one(T)
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
by `-theta` to align with the principal axes.

The model is evaluated as 

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
peak(model::GaussianPRF{T}) where T =
    model.flux * erf(sqrt(T(log(2))) / model.x_fwhm) * erf(sqrt(T(log(2))) / model.y_fwhm) + model.bkg
effective_area(model::GaussianPRF{T}) where T = π * model.x_fwhm * model.y_fwhm / T(2 * log(2))

function evaluate(model::GaussianPRF{T}, px, py) where T
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
function evaluate_fg(model::GaussianPRF{T}, px, py) where T
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
    df_dx    = fl4 * (cs * dEx_du * Ey - sn * dEy_dv * Ex)
    # ∂f/∂y: ∂u/∂y = -sn, ∂v/∂y = -cs
    #   ∂f/∂y = fl4*((-dEx_du)*(-sn)*Ey + Ex*(-dEy_dv)*(-cs)) = fl4*(dEx_du*sn*Ey + dEy_dv*cs*Ex)
    df_dy    = fl4 * (sn * dEx_du * Ey + cs * dEy_dv * Ex)
    # ∂f/∂x_fwhm: same form but with u replacing dx
    df_dax   = fl4 / model.x_fwhm * (Gxm * u_m - Gxp * u_p) * Ey
    df_day   = fl4 / model.y_fwhm * Ex * (Gym * v_m - Gyp * v_p)
    # ∂f/∂θ: ∂u/∂θ = d*v, ∂v/∂θ = -d*u, where d = deg2rad(1) = π / 180
    #   ∂f/∂θ = fl4*((-dEx_du)*d*v*Ey + Ex*(-dEy_dv)*(-d*u)) = fl4*d*(dEy_dv*u*Ex - dEx_du*v*Ey)
    df_dtheta = fl4 * deg2rad(one(T)) * (dEy_dv * u * Ex - dEx_du * v * Ey)
    df_dflux = Ex * Ey / 4
    df_dbkg  = one(T)
    G = SA[df_dx, df_dy, df_dax, df_day, df_dtheta, df_dflux, df_dbkg]
    return f, G
end
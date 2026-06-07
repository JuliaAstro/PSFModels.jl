function _as_oversampling(oversampling)
    # Accept scalar oversampling as the same factor along both axes.
    if oversampling isa Integer
        oversampling > 0 || throw(ArgumentError("`oversampling` must be positive"))
        return (Int(oversampling), Int(oversampling))
    elseif oversampling isa Tuple || oversampling isa AbstractVector
        # Accept explicit axis factors for anisotropic sampling.
        length(oversampling) == 2 ||
            throw(ArgumentError("`oversampling` must be an integer or a length-2 tuple/vector of integers"))
        all(x -> x isa Integer, oversampling) ||
            throw(ArgumentError("length-2 `oversampling` values must be integers"))
        sx, sy = Int(oversampling[1]), Int(oversampling[2])
        sx > 0 && sy > 0 || throw(ArgumentError("`oversampling` values must be positive"))
        return (sx, sy)
    else
        throw(ArgumentError("`oversampling` must be an integer or a length-2 tuple/vector of integers"))
    end
end

"""
    ImagePSF(data; x=0, y=0, flux=1, bkg=0, origin=nothing,
             oversampling=1, fill_value=0, normalize=false)

Empirical image-backed PSF/ePSF model. `data` is a finite 2D sampled PSF image;
`x` and `y` are the source centroid in detector-pixel coordinates, `flux` is a
multiplicative scale, and `bkg` is an additive scalar background. `oversampling`
is the integer sampling factor of `data` relative to detector pixels.

The model evaluates the tabulated PSF with separable 4x4 bicubic interpolation.
Values outside the tabulated image are `flux * fill_value + bkg`.

# Notes
## Centroid and origin
`x` and `y` give the center of the star *in the image being modeled*. `origin` locates the
stellar centroid inside the tabulated PSF image, in 1-based array coordinates.
When evaluating at detector-pixel coordinates `(px, py)` for a star whose putative
center is at `(x, y)`, the tabulated PSF is sampled at

```julia
u = oversampling[1] * (px - x) + origin[1]
v = oversampling[2] * (py - y) + origin[2]
```

Thus changing `x` and `y` moves the source through the detector image, while
changing `origin` changes which point in `data` is interpreted as the PSF
center. For ordinary fitted sources, `origin` should remain fixed and only
`x`, `y`, `flux`, and optionally `bkg` should vary.

## Oversampling
Oversampling is useful when many stars sample different subpixel phases of the
same effective PSF. Higher oversampling requires more well-centered,
isolated stars; otherwise many oversampled cells will be
weakly constrained and must be filled by interpolation. A few dozen stars can be
enough for modest factors such as 2x in clean simulations, but real undersampled
data generally benefit from hundreds of good stars when available. Very high
oversampling with too few stars can produce a smooth-looking but poorly
constrained model.

## Normalization
If `normalize=true`, the input `data` is rescaled so that `sum(data) ==
prod(oversampling)`. This is the convention used by the evaluator for a
unit-flux ePSF: a fitted `flux` remains on the detector-pixel flux scale. Use
`normalize=true` when constructing an `ImagePSF` from an arbitrary image stamp
or stack whose sum is not already in this convention. Use `normalize=false`
when `data` is already normalized, when preserving the exact tabulated values is
important, or when intentionally using `flux` as a pure amplitude multiplier
rather than a physical total flux.
"""
struct ImagePSF{T, A <: AbstractMatrix{T}} <: AbstractPSFModel{T}
    data::A
    x::T
    y::T
    flux::T
    bkg::T
    origin::Tuple{T, T}
    oversampling::Tuple{Int, Int}
    fill_value::T
end

function ImagePSF(
        data::AbstractMatrix;
        x = 0,
        y = 0,
        flux = 1,
        bkg = 0,
        origin = nothing,
        oversampling = 1,
        fill_value = 0,
        normalize::Bool = false
    )
    # Validate the tabulated PSF before building interpolation metadata.
    ndims(data) == 2 || throw(ArgumentError("`data` must be a 2D matrix"))
    all(≥(4), size(data)) || throw(ArgumentError("each axis of `data` must have at least 4 samples"))
    all(isfinite, data) || throw(ArgumentError("all elements of `data` must be finite"))
    os = _as_oversampling(oversampling)

    # Promote model parameters and grid values to one floating-point type.
    T = promote_type(eltype(data), typeof(x), typeof(y), typeof(flux), typeof(bkg), typeof(fill_value))
    if !isnothing(origin)
        length(origin) == 2 || throw(ArgumentError("`origin` must have two elements"))
        T = promote_type(T, typeof(origin[1]), typeof(origin[2]))
    end
    T = float(T)

    # Allocating new matrix here is not *strictly* necessary but it is safe;
    # external modifications to the input `data` after model
    # construction would be surprising and could cause hard-to-debug issues.
    dataT = Matrix{T}(data)
    # Optionally normalize so the oversampled grid sums to the sampling area.
    if normalize
        s = sum(dataT)
        isfinite(s) && s > zero(T) || throw(ArgumentError("cannot normalize PSF data with non-positive sum"))
        dataT .*= T(os[1] * os[2]) / s
    end

    # Default origin is the central grid sample; otherwise use the caller's origin.
    org = if isnothing(origin)
        (T((size(dataT, 1) + 1) / 2), T((size(dataT, 2) + 1) / 2))
    else
        ox, oy = T(origin[1]), T(origin[2])
        isfinite(ox) && isfinite(oy) || throw(ArgumentError("`origin` must be finite"))
        (ox, oy)
    end

    return ImagePSF{T, typeof(dataT)}(dataT, T(x), T(y), T(flux), T(bkg), org, os, T(fill_value))
end

ImagePSF(data::AbstractMatrix, x, y, flux, bkg; kwargs...) =
    ImagePSF(data; x, y, flux, bkg, kwargs...)

ConstructionBase.getproperties(model::ImagePSF) = (x = model.x, y = model.y, flux = model.flux, bkg = model.bkg)

function ConstructionBase.setproperties(model::ImagePSF{T, S}, patch::NamedTuple) where {T, S}
    # Only fit parameters can change; PSF data stay fixed.
    x = haskey(patch, :x) ? T(patch.x) : model.x
    y = haskey(patch, :y) ? T(patch.y) : model.y
    flux = haskey(patch, :flux) ? T(patch.flux) : model.flux
    bkg = haskey(patch, :bkg) ? T(patch.bkg) : model.bkg
    return ImagePSF{T, S}(
        model.data, x, y, flux, bkg, model.origin,
        model.oversampling, model.fill_value
    )
end

theta(model::ImagePSF{T}) where {T} = zero(T)

function extent(model::ImagePSF{T}) where {T}
    # Convert tabulated grid bounds back into detector-pixel coordinates.
    nx, ny = size(model.data)
    sx, sy = model.oversampling
    ox, oy = model.origin
    return (
        (model.x - (ox - one(T)) / sx, model.x + (T(nx) - ox) / sx),
        (model.y - (oy - one(T)) / sy, model.y + (T(ny) - oy) / sy),
    )
end
extent(model::ImagePSF, _) = extent(model) # ignore extra argument (fwhm_factor) if provided.

function effective_area(model::ImagePSF{T}) where {T}
    # Approximate effective area from the normalized discrete ePSF samples.
    denom = sum(abs2, model.data)
    return denom == zero(T) ? T(Inf) : T(model.oversampling[1] * model.oversampling[2]) / denom
end

@inline function _cubic4(f1, f2, f3, f4, t)
    # Build cubic-convolution coefficients from four adjacent samples.
    c1 = (f3 - f1) / 2
    c4 = f3 - f2 - c1
    c2 = 3 * c4 - (f4 - f2) / 2 + c1
    c3 = c4 - c2
    # Evaluate both the interpolated value and its derivative at t.
    c4t = t * c3
    value = t * (t * (c4t + c2) + c1) + f2
    deriv = t * (3 * c4t + 2 * c2) + c1
    return value, deriv
end

@doc raw"""
    bicubic_interpolate(data, x, y; fill_value=0)

Evaluate `data` at 1-based fractional array coordinates `(x, y)` with a
separable 4x4 cubic-convolution interpolant. Returns `(value, dfdx, dfdy)`,
where `dfdx` and `dfdy` are derivatives with respect to the input grid
coordinates.

For a coordinate inside the array, let `lx = floor(x)`, `ly = floor(y)`,
`dx = x - lx`, and `dy = y - ly`. The interpolation point is treated as lying
between grid columns `lx` and `lx + 1`, and between rows `ly` and `ly + 1`.
For each of the four neighboring rows `ly-1:ly+2`, the four samples
`data[lx-1:lx+2, row]` define a cubic polynomial in `dx`:

```math
p(t) = a_0 + a_1 t + a_2 t^2 + a_3 t^3,\qquad 0 \le t \le 1.
```

The coefficients are chosen by the cubic-convolution constraints

```math
p(0) = f_2,\quad p(1) = f_3,\quad
p'(0) = \frac{f_3 - f_1}{2},\quad
p'(1) = \frac{f_4 - f_2}{2},
```

where `(f_1, f_2, f_3, f_4)` are the four row samples. This gives one
interpolated value and one x-derivative for each row. The four row values are
then interpolated with the same cubic-convolution rule in `dy` to produce the
final value and `dfdy`. The four row-wise x-derivatives are also interpolated in
`dy` to produce `dfdx`.

Boundary conditions and assumptions:
- `data` must have at least four samples along each axis.
- Coordinates outside the closed array bounds return `(fill_value, 0, 0)`.
- Inside the array, the 4x4 stencil is clamped at the nearest array edge. This
  is equivalent to constant extrapolation of edge samples for the stencil only.
- The interpolant is intended for regularly spaced grid samples. `x` and `y`
  are array coordinates, not detector-pixel coordinates.
"""
function bicubic_interpolate(data::AbstractMatrix, x, y; fill_value = zero(eltype(data)))
    # Work in a promoted floating type and reject out-of-grid samples early.
    T = promote_type(eltype(data), typeof(x), typeof(y), typeof(fill_value))
    T = float(T)
    nx, ny = size(data)
    nx ≥ 4 && ny ≥ 4 || throw(ArgumentError("`data` must have at least four samples along each axis"))
    xx, yy = T(x), T(y)
    if !(isfinite(xx) && isfinite(yy)) || xx < one(T) || xx > T(nx) || yy < one(T) || yy > T(ny)
        return T(fill_value), zero(T), zero(T)
    end

    # Locate the lower-left grid cell and fractional offset inside it.
    lx = clamp(floor(Int, xx), 1, nx - 1)
    ly = clamp(floor(Int, yy), 1, ny - 1)
    dx = xx - T(lx)
    dy = yy - T(ly)

    # Clamp the 4x4 stencil once so edge samples are reused directly below.
    ix1 = clamp(lx - 1, 1, nx)
    ix2 = lx
    ix3 = lx + 1
    ix4 = clamp(lx + 2, 1, nx)
    iy1 = clamp(ly - 1, 1, ny)
    iy2 = ly
    iy3 = ly + 1
    iy4 = clamp(ly + 2, 1, ny)

    # Interpolate each row in x and keep the row-wise x derivatives.
    @inbounds begin
        row1, drow1dx = _cubic4(
            T(data[ix1, iy1]), T(data[ix2, iy1]),
            T(data[ix3, iy1]), T(data[ix4, iy1]), dx
        )
        row2, drow2dx = _cubic4(
            T(data[ix1, iy2]), T(data[ix2, iy2]),
            T(data[ix3, iy2]), T(data[ix4, iy2]), dx
        )
        row3, drow3dx = _cubic4(
            T(data[ix1, iy3]), T(data[ix2, iy3]),
            T(data[ix3, iy3]), T(data[ix4, iy3]), dx
        )
        row4, drow4dx = _cubic4(
            T(data[ix1, iy4]), T(data[ix2, iy4]),
            T(data[ix3, iy4]), T(data[ix4, iy4]), dx
        )
    end

    # Interpolate those row values in y, including both first derivatives.
    value, dfdy = _cubic4(row1, row2, row3, row4, dy)
    dfdx, _ = _cubic4(drow1dx, drow2dx, drow3dx, drow4dx, dy)
    return value, dfdx, dfdy
end

function evaluate(model::ImagePSF{T}, px, py) where {T}
    # Map detector-pixel coordinates into the oversampled PSF grid.
    sx, sy = model.oversampling
    u = T(sx) * (T(px) - model.x) + model.origin[1]
    v = T(sy) * (T(py) - model.y) + model.origin[2]
    # Scale the unit-flux ePSF sample and add scalar background.
    p, _, _ = bicubic_interpolate(model.data, u, v; fill_value = model.fill_value)
    return muladd(model.flux, T(p), model.bkg)
end

function evaluate_fg(model::ImagePSF{T}, px, py) where {T}
    # Evaluate the ePSF and its grid-space derivatives at this detector pixel.
    sx, sy = model.oversampling
    u = T(sx) * (T(px) - model.x) + model.origin[1]
    v = T(sy) * (T(py) - model.y) + model.origin[2]
    p, dpdu, dpdv = bicubic_interpolate(model.data, u, v; fill_value = model.fill_value)
    # Apply the chain rule for centroid, flux, and background parameters.
    profile = T(p)
    f = muladd(model.flux, profile, model.bkg)
    df_dx = -model.flux * T(sx) * T(dpdu)
    df_dy = -model.flux * T(sy) * T(dpdv)
    df_dflux = profile
    df_dbkg = one(T)
    return f, SA[df_dx, df_dy, df_dflux, df_dbkg]
end

"""
    ImagePSFBuildResult

Metadata returned by `fit(ImagePSF, ...)`.
"""
struct ImagePSFBuildResult{T}
    psf::ImagePSF{T}
    x::Vector{T}
    y::Vector{T}
    flux::Vector{T}
    bkg::Vector{T}
    used::BitVector
    converged::BitVector
    iterations::Int
    costs::Vector{T}
end

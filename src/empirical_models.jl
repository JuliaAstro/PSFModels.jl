"""The quartic smoothing kernel of [`Anderson2000`](@cite) (Eq. 8) used by default when `smooth=true` in `fit(ImagePSF, ...)`."""
const _QUARTIC_SMOOTHING_KERNEL = (
    (0.041632, -0.080816, 0.078368, -0.080816, 0.041632),
    (-0.080816, -0.019592, 0.200816, -0.019592, -0.080816),
    (0.078368, 0.200816, 0.441632, 0.200816, 0.078368),
    (-0.080816, -0.019592, 0.200816, -0.019592, -0.080816),
    (0.041632, -0.080816, 0.078368, -0.080816, 0.041632),
)

function _as_oversampling(oversampling)
    # Accept scalar oversampling as the same factor along both axes.
    if oversampling isa Integer
        oversampling > 0 || throw(ArgumentError("`oversampling` must be positive"))
        return (Int(oversampling), Int(oversampling))
    # Accept explicit axis factors for anisotropic sampling.
    elseif oversampling isa Tuple || oversampling isa AbstractVector
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
        model.data,
        x,
        y,
        flux,
        bkg,
        model.origin,
        model.oversampling,
        model.fill_value
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
        row1, drow1dx = _cubic4(T(data[ix1, iy1]), T(data[ix2, iy1]), 
                                T(data[ix3, iy1]), T(data[ix4, iy1]), dx)
        row2, drow2dx = _cubic4(T(data[ix1, iy2]), T(data[ix2, iy2]), 
                                T(data[ix3, iy2]), T(data[ix4, iy2]), dx)
        row3, drow3dx = _cubic4(T(data[ix1, iy3]), T(data[ix2, iy3]), 
                                T(data[ix3, iy3]), T(data[ix4, iy3]), dx)
        row4, drow4dx = _cubic4(T(data[ix1, iy4]), T(data[ix2, iy4]), 
                                T(data[ix3, iy4]), T(data[ix4, iy4]), dx)
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

mutable struct _EmpiricalStar{T}
    inds::Tuple{UnitRange{Int}, UnitRange{Int}}
    pixels::Vector{CartesianIndex{2}}
    x::T
    y::T
    flux::T
    bkg::T
    used::Bool
    converged::Bool
    cost::T
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

@inline function _round_half_away(x)
    return x ≥ zero(x) ? floor(Int, x + oftype(x, 0.5)) : ceil(Int, x - oftype(x, 0.5))
end

@inline function _pixel_wholly_inside(i, j, x, y, r)
    return (abs(i - x) + 0.5)^2 + (abs(j - y) + 0.5)^2 ≤ r^2
end

function _edge_background(image, inds)
    # Use finite edge pixels as a simple local sky estimate for the cutout.
    xs, ys = inds
    vals = Float64[]
    for i in xs, j in ys
        if i == first(xs) || i == last(xs) || j == first(ys) || j == last(ys)
            v = image[i, j]
            isfinite(v) && push!(vals, float(v))
        end
    end
    return isempty(vals) ? zero(float(eltype(image))) : median(vals)
end

function _initialize_star!(star::_EmpiricalStar, image)
    # Estimate local sky first so flux uses background-subtracted pixels.
    bkg = _edge_background(image, star.inds)
    vals = typeof(star.x)[]
    for idx in star.pixels
        v = image[idx]
        isfinite(v) && push!(vals, max(zero(typeof(star.x)), typeof(star.x)(v - bkg)))
    end
    # Winsorize the brightest pixels so hot defects do not set the aperture flux.
    sort!(vals)
    cap = isempty(vals) ? zero(typeof(star.x)) : vals[clamp(ceil(Int, 0.98 * length(vals)), 1, length(vals))]
    flux = sum(v -> min(v, cap), vals; init = zero(typeof(star.x)))
    # Fall back to a peak estimate if the aperture sum is unusable.
    if !(isfinite(flux) && flux > eps(typeof(flux)))
        flux = max(eps(typeof(flux)), typeof(flux)(maximum(image[star.inds...] .- bkg)))
    end
    star.bkg = typeof(star.bkg)(bkg)
    star.flux = typeof(star.flux)(flux)
    return star
end

function _stars_from_xy(image, x, y, fit_rad; drop_edge::Bool)
    # Normalize coordinate inputs and choose a working floating type.
    xs = collect(x)
    ys = collect(y)
    length(xs) == length(ys) || throw(ArgumentError("`x` and `y` must have the same length"))
    T = promote_type(eltype(image), eltype(xs), eltype(ys), typeof(fit_rad))
    T = T <: Integer ? Float64 : float(T)
    ax, ay = axes(image)
    stars = _EmpiricalStar{T}[]
    for k in eachindex(xs, ys)
        # Convert each initial center into a square candidate cutout.
        xk, yk = T(xs[k]), T(ys[k])
        xlo = floor(Int, xk - fit_rad)
        xhi = ceil(Int, xk + fit_rad)
        ylo = floor(Int, yk - fit_rad)
        yhi = ceil(Int, yk + fit_rad)
        inside = first(ax) ≤ xlo && xhi ≤ last(ax) && first(ay) ≤ ylo && yhi ≤ last(ay)
        if !inside
            drop_edge && continue
            throw(ArgumentError("fit cutout for star $k extends outside the image"))
        end

        # Keep only pixels whose full pixel square lies inside the fit radius.
        pix = CartesianIndex{2}[]
        for i in xlo:xhi, j in ylo:yhi
            _pixel_wholly_inside(i, j, xk, yk, fit_rad) && push!(pix, CartesianIndex(i, j))
        end
        isempty(pix) && throw(ArgumentError("`fit_rad` leaves no usable pixels for star $k"))
        # Initialize per-star sky and flux estimates before stacking.
        star = _EmpiricalStar((xlo:xhi, ylo:yhi), pix, xk, yk, one(T), zero(T), true, false, T(Inf))
        _initialize_star!(star, image)
        push!(stars, star)
    end
    isempty(stars) && throw(ArgumentError("no usable stars remain after cutout extraction"))
    return stars
end

function _normalize_cutout_inds(inds)
    # Treat one `(xrange, yrange)` tuple as a single cutout, not two cutouts.
    if inds isa Tuple && length(inds) == 2 && inds[1] isa AbstractUnitRange && inds[2] isa AbstractUnitRange
        return (inds,)
    end
    return Tuple(inds)
end

function _stars_from_inds(image, inds; x = nothing, y = nothing)
    # Normalize cutout input and validate optional initial centers.
    cutouts = _normalize_cutout_inds(inds)
    n = length(cutouts)
    if !isnothing(x) || !isnothing(y)
        (isnothing(x) || isnothing(y)) && throw(ArgumentError("`x` and `y` must be provided together"))
        length(x) == n && length(y) == n || throw(ArgumentError("`x` and `y` must match `inds` length"))
    end
    T = promote_type(eltype(image), isnothing(x) ? Float64 : eltype(x), isnothing(y) ? Float64 : eltype(y))
    T = T <: Integer ? Float64 : float(T)
    stars = _EmpiricalStar{T}[]
    for k in 1:n
        # Convert each cutout to concrete integer ranges within the image.
        c = cutouts[k]
        length(c) == 2 || throw(ArgumentError("each cutout must be a tuple of two ranges"))
        xr = Int(first(c[1])):Int(last(c[1]))
        yr = Int(first(c[2])):Int(last(c[2]))
        first(axes(image, 1)) ≤ first(xr) && last(xr) ≤ last(axes(image, 1)) ||
            throw(ArgumentError("x range for cutout $k is outside the image"))
        first(axes(image, 2)) ≤ first(yr) && last(yr) ≤ last(axes(image, 2)) ||
            throw(ArgumentError("y range for cutout $k is outside the image"))
        # Use caller-supplied centers or the geometric midpoint of the cutout.
        xk = isnothing(x) ? T((first(xr) + last(xr)) / 2) : T(x[k])
        yk = isnothing(y) ? T((first(yr) + last(yr)) / 2) : T(y[k])
        pix = vec(collect(CartesianIndices((xr, yr))))
        # Initialize each cutout as a fitting star.
        star = _EmpiricalStar((xr, yr), pix, xk, yk, one(T), zero(T), true, false, T(Inf))
        _initialize_star!(star, image)
        push!(stars, star)
    end
    isempty(stars) && throw(ArgumentError("`inds` must contain at least one cutout"))
    return stars
end

function _robust_location(vals::Vector{T}, sigma_clip) where {T}
    # Median is the baseline stack statistic.
    isempty(vals) && return T(NaN)
    med = median(vals)
    if length(vals) < 4 || isnothing(sigma_clip)
        return med
    end
    # Sigma-clip around the median when there are enough samples.
    σ = T(1.4826) * median(abs.(vals .- med))
    if !(isfinite(σ) && σ > eps(T))
        return med
    end
    keep = T[]
    c = T(sigma_clip) * σ
    for v in vals
        abs(v - med) ≤ c && push!(keep, v)
    end
    return isempty(keep) ? med : median(keep)
end

function _nearest4_finite_indices(data, i, j; along_x::Bool)
    # Search outward from the requested index until four finite samples exist.
    n = along_x ? size(data, 1) : size(data, 2)
    center = along_x ? i : j
    candidates = Int[]
    for r in 0:n
        for sgn in (-1, 1)
            idx = center + sgn * r
            1 ≤ idx ≤ n || continue
            if r == 0 && sgn == 1
                continue
            end
            v = along_x ? data[idx, j] : data[i, idx]
            if isfinite(v)
                push!(candidates, idx)
                length(candidates) == 4 && return candidates
            end
        end
    end
    return candidates
end

function _has_finite_bracket(data, i, j; along_x::Bool)
    # Require finite samples on both sides so infill interpolates rather than
    # extrapolates beyond the sampled stellar footprint.
    n = along_x ? size(data, 1) : size(data, 2)
    center = along_x ? i : j
    has_lo = false
    has_hi = false
    for idx in 1:n
        idx == center && continue
        v = along_x ? data[idx, j] : data[i, idx]
        isfinite(v) || continue
        has_lo |= idx < center
        has_hi |= idx > center
        has_lo && has_hi && return true
    end
    return false
end

function _cubic_lagrange(xs, ys, x0)
    # Evaluate the cubic Lagrange interpolant through four arbitrary samples.
    # For infill, nearby finite samples may not form a complete adjacent stencil,
    # so this is more useful than the regular-grid bicubic evaluator.
    T = promote_type(eltype(xs), eltype(ys), typeof(x0))
    out = zero(T)
    for k in 1:4
        w = one(T)
        xk = T(xs[k])
        for m in 1:4
            m == k && continue
            w *= (T(x0) - T(xs[m])) / (xk - T(xs[m]))
        end
        out += T(ys[k]) * w
    end
    return out
end

"""
    _fill_missing_bicubic!(data; maxiter::Int = 6))

Fill non-finite holes in an empirical ePSF stack in place. Interior holes are
filled first with bracketed separable cubic interpolation, unsupported border
holes then use a local exponential radial extrapolation 
(see 4.2.1 of [Anderson2000](@citet)) and any remaining holes fall back to a
local median if possible, and zero fill otherwise. `_stack_psf` uses
this after robustly combining stellar samples into the *initial* oversampled grid
and before smoothing, recentering, and normalization. `maxiter` controls the
maximum number of passes for the bicubic and median stages; six is a reasonable
default to allow infill growth without risking infinite loops.
"""
function _fill_missing_bicubic!(data; maxiter::Int = 6)
    # First try a separable cubic fill using nearby finite rows and columns.
    T = eltype(data)
    nx, ny = size(data)
    # Skip all infill work when the stack already has no holes.
    all(isfinite, data) && return data

    fillable = falses(nx, ny)
    for j in 1:ny, i in 1:nx
        # Mark only holes with finite support on both sides along both axes.
        if !isfinite(data[i, j])
            fillable[i, j] = _has_finite_bracket(data, i, j; along_x = true) &&
                _has_finite_bracket(data, i, j; along_x = false)
        end
    end
    if count(fillable) / (nx * ny) > 0.1
        @warn "more than 10% of the ePSF grid is missing; infill may be unreliable"
    end
    # Iterate because filling one interior hole can provide the finite support
    # needed to bicubically fill adjacent holes; six passes bounds the growth.
    for _ in 1:maxiter
        changed = false
        old = copy(data)
        for j in 1:ny, i in 1:nx
            isfinite(old[i, j]) && continue
            fillable[i, j] || continue
            yinds = _nearest4_finite_indices(old, i, j; along_x = false)
            length(yinds) == 4 || continue
            any(yy -> yy < j, yinds) && any(yy -> yy > j, yinds) || continue
            row_vals = T[]
            good_y = Int[]
            for yy in yinds
                _has_finite_bracket(old, i, yy; along_x = true) || continue
                xinds = _nearest4_finite_indices(old, i, yy; along_x = true)
                length(xinds) == 4 || continue
                any(xx -> xx < i, xinds) && any(xx -> xx > i, xinds) || continue
                vals = T[old[xx, yy] for xx in xinds]
                push!(row_vals, _cubic_lagrange(xinds, vals, i))
                push!(good_y, yy)
            end
            if length(row_vals) == 4
                data[i, j] = _cubic_lagrange(good_y, row_vals, j)
                changed = true
            end
        end
        changed || break
    end

    # Extrapolate border holes with Anderson's local exponential radial tail.
    cx = T((nx + 1) / 2)
    cy = T((ny + 1) / 2)
    old = copy(data)
    for j in 1:ny, i in 1:nx
        isfinite(old[i, j]) && continue
        fillable[i, j] && continue
        rs = T[]
        logs = T[]
        vals = T[]
        for jj in max(1, j - 2):min(ny, j + 2), ii in max(1, i - 2):min(nx, i + 2)
            # Use only positive finite samples because the fit is done in log space.
            v = old[ii, jj]
            if isfinite(v) && v > zero(T)
                push!(rs, hypot(T(ii) - cx, T(jj) - cy))
                push!(logs, log(T(v)))
                push!(vals, T(v))
            end
        end
        length(rs) ≥ 2 || continue

        # Fit log(value) = a + b * radius using the local positive support.
        rmean = sum(rs) / length(rs)
        lmean = sum(logs) / length(logs)
        denom = sum(abs2, r - rmean for r in rs; init = zero(T))
        denom > eps(T) || continue
        slope = sum(
            (rs[k] - rmean) * (logs[k] - lmean) for k in eachindex(rs);
            init = zero(T)
        ) / denom
        slope < zero(T) || continue
        intercept = lmean - slope * rmean

        # Accept only outward, decaying extrapolations to avoid false border peaks.
        target_r = hypot(T(i) - cx, T(j) - cy)
        target_r ≥ minimum(rs) || continue
        value = exp(intercept + slope * target_r)
        if isfinite(value) && zero(T) < value ≤ maximum(vals)
            data[i, j] = value
        end
    end

    # Fill remaining holes with local neighbor medians.
    for _ in 1:maxiter
        changed = false
        old = copy(data)
        for j in 1:ny, i in 1:nx
            isfinite(old[i, j]) && continue
            vals = T[]
            for jj in max(1, j - 2):min(ny, j + 2), ii in max(1, i - 2):min(nx, i + 2)
                isfinite(old[ii, jj]) && push!(vals, old[ii, jj])
            end
            if !isempty(vals)
                data[i, j] = median(vals)
                changed = true
            end
        end
        changed || break
    end
    # Last resort: make sure the PSF grid is fully finite.
    for idx in eachindex(data)
        isfinite(data[idx]) || (data[idx] = zero(T))
    end
    return data
end

function _smooth_quartic(data::AbstractMatrix{T}) where {T}
    # Apply the fixed quartic ePSF smoothing kernel with clamped edges.
    # Manual 2-D discrete convolution with fixed 5x5 kernel.
    nx, ny = size(data)
    out = similar(data)
    for j in 1:ny, i in 1:nx
        acc = zero(T)
        for kj in 1:5, ki in 1:5
            ii = clamp(i + ki - 3, 1, nx)
            jj = clamp(j + kj - 3, 1, ny)
            acc += T(_QUARTIC_SMOOTHING_KERNEL[kj][ki]) * data[ii, jj]
        end
        out[i, j] = acc
    end
    return out
end

function _recenter_data(data, origin, oversampling)
    # Use positive ePSF samples to estimate the current grid centroid.
    T = eltype(data)
    weights = max.(data, zero(T))
    s = sum(weights)
    s > eps(T) || return data
    cx = zero(T)
    cy = zero(T)
    for j in axes(data, 2), i in axes(data, 1)
        w = weights[i, j]
        cx += T(i) * w
        cy += T(j) * w
    end
    cx /= s
    cy /= s
    dx = cx - T(origin[1])
    dy = cy - T(origin[2])
    if abs(dx) < T(1.0e-6) && abs(dy) < T(1.0e-6)
        return data
    end
    # Ignore large shifts; those usually signal a bad stack, not a recentering need.
    if abs(dx) > T(2 * oversampling[1]) || abs(dy) > T(2 * oversampling[2])
        return data
    end
    # Shift the sampled grid back onto the requested origin.
    out = similar(data)
    for j in axes(data, 2), i in axes(data, 1)
        out[i, j], _, _ = bicubic_interpolate(data, T(i) + dx, T(j) + dy; fill_value = zero(T))
    end
    return out
end

function _normalize_psf_data!(data, oversampling; clip_negative::Bool)
    # Optionally remove negative stack artifacts before normalization.
    T = eltype(data)
    if clip_negative
        for idx in eachindex(data)
            data[idx] = max(zero(T), data[idx])
        end
    end
    # Enforce the ImagePSF convention: sum equals oversampling area.
    s = sum(data)
    if !(isfinite(s) && s > eps(T))
        throw(ArgumentError("empirical PSF stack has non-positive finite sum"))
    end
    data .*= T(oversampling[1] * oversampling[2]) / s
    return data
end

function _stack_psf(
        image,
        stars,
        shape,
        origin,
        oversampling;
        sigma_clip,
        sample_clip,
        smooth,
        recenter::Bool,
        clip_negative::Bool,
        badmask = nothing
    )
    # Allocate one sample list per oversampled ePSF grid cell.
    T = promote_type(float(eltype(image)), typeof(stars[1].x))
    nx, ny = shape
    cells = [T[] for _ in 1:(nx * ny)]
    sx, sy = oversampling
    ox, oy = origin
    for star in stars
        star.used || continue
        isfinite(star.flux) && star.flux > eps(T) || continue
        # Project each valid, normalized stellar pixel into ePSF-grid space.
        for idx in star.pixels
            isnothing(badmask) || !badmask[idx] || continue
            v = image[idx]
            isfinite(v) || continue
            gx = _round_half_away(T(ox) + T(sx) * (T(idx[1]) - star.x))
            gy = _round_half_away(T(oy) + T(sy) * (T(idx[2]) - star.y))
            1 ≤ gx ≤ nx && 1 ≤ gy ≤ ny || continue
            # Drop impossible normalized samples from severe pixel defects.
            sample = T(v - star.bkg) / star.flux
            if !isnothing(sample_clip) && abs(sample) > T(sample_clip)
                continue
            end
            push!(cells[gx + (gy - 1) * nx], sample)
        end
    end

    # Robustly combine samples in each grid cell.
    data = fill(T(NaN), nx, ny)
    for gy in 1:ny, gx in 1:nx
        vals = cells[gx + (gy - 1) * nx]
        isempty(vals) || (data[gx, gy] = _robust_location(vals, sigma_clip))
    end
    # Fill holes, optionally smooth, optionally recenter, then normalize.
    _fill_missing_bicubic!(data)
    if smooth === true || smooth === :quartic
        data = _smooth_quartic(data)
    elseif smooth === false || isnothing(smooth)
        nothing
    else
        throw(ArgumentError("`smooth` must be true, false, `:quartic`, or `nothing`"))
    end
    recenter && (data = _recenter_data(data, origin, oversampling))
    _normalize_psf_data!(data, oversampling; clip_negative)
    return data
end

function _valid_star_pixels(star, image, badmask)
    # Return finite unmasked pixels for a per-star LM fit.
    pix = CartesianIndex{2}[]
    for idx in star.pixels
        isnothing(badmask) || !badmask[idx] || continue
        isfinite(image[idx]) && push!(pix, idx)
    end
    return pix
end

function _fit_star!(
        star::_EmpiricalStar{T},
        psf::ImagePSF,
        image;
        badmask = nothing,
        reweight = TukeyLoss(),
        star_max_iter::Integer = 100,
        fit_bkg::Bool = false,
        show_trace::Bool = false,
        kwargs...
    ) where {T}
    # Select finite, unmasked observations and make sure the fit is constrained.
    pix = _valid_star_pixels(star, image, badmask)
    nparams = fit_bkg ? 4 : 3
    if length(pix) <= nparams
        star.used = false
        star.converged = false
        return nothing
    end
    x0 = fit_bkg ? T[star.x, star.y, star.flux, star.bkg] : T[star.x, star.y, star.flux]

    # Stream residuals and Jacobians into LM normal equations.
    function accum!(A::AbstractMatrix{FT}, b::AbstractVector{FT}, residuals::AbstractVector{FT}, x::AbstractVector{FT}, weights) where {FT}
        fill!(A, zero(FT))
        fill!(b, zero(FT))
        cost = zero(FT)
        # Rebuild the ImagePSF with trial source parameters.
        model = ImagePSF(
            psf.data;
            x = x[1],
            y = x[2],
            flux = x[3],
            bkg = fit_bkg ? x[4] : star.bkg,
            origin = psf.origin,
            oversampling = psf.oversampling,
            fill_value = psf.fill_value
        )
        # Accumulate weighted least-squares terms pixel by pixel.
        @inbounds for k in eachindex(pix)
            idx = pix[k]
            w = isnothing(weights) ? one(FT) : FT(weights[k])
            fval, g = evaluate_fg(model, idx)
            r = FT(fval) - FT(image[idx])
            residuals[k] = r
            wr = w * r
            cost = muladd(wr, r, cost)
            for j in 1:nparams
                gj = FT(g[j])
                b[j] = muladd(wr, gj, b[j])
                for i in 1:nparams
                    A[i, j] = muladd(w * FT(g[i]), gj, A[i, j])
                end
            end
        end
        return cost
    end

    # Run LM/IRLS on this star's centroid and flux parameters.
    problem = LMProblem(x0, length(pix), accum!)
    result = lm_irls(
        problem;
        max_iter = star_max_iter,
        reweight,
        show_trace,
        kwargs...
    )
    # Accept only finite positive-flux fits; otherwise disable this star.
    p = result.minimizer
    if all(isfinite, p) && p[3] > eps(T)
        star.x = T(p[1])
        star.y = T(p[2])
        star.flux = T(p[3])
        fit_bkg && (star.bkg = T(p[4]))
        star.converged = result.converged
        star.cost = T(result.minimum)
        star.used = result.converged
    else
        star.used = false
        star.converged = false
        star.cost = T(Inf)
    end
    return result
end

function _fit_imagepsf(
        image::AbstractMatrix,
        stars::Vector{_EmpiricalStar{T}};
        oversampling = 1,
        psf_radius = nothing,
        maxiters::Integer = 6,
        centroid_tol::Real = 1.0e-3,
        sigma_clip = 4.0,
        sample_clip = 2.0,
        smooth = true,
        recenter::Bool = true,
        anchor_centroids::Bool = true,
        clip_negative::Bool = true,
        fill_value = 0,
        badmask = nothing,
        fit_bkg::Bool = false,
        reweight = TukeyLoss(),
        star_max_iter::Integer = 100,
        show_trace::Bool = false,
        kwargs...
    ) where {T}
    # Validate optional masks and normalize oversampling metadata.
    if !isnothing(badmask)
        size(badmask) == size(image) || throw(ArgumentError("`badmask` must have the same size as `image`"))
    end
    os = _as_oversampling(oversampling)
    # Choose an odd supersampled PSF grid with a central origin.
    radius = if isnothing(psf_radius)
        maximum(stars) do star
            maximum(max(abs(idx[1] - star.x), abs(idx[2] - star.y)) for idx in star.pixels)
        end
    else
        psf_radius
    end
    nx = 2 * ceil(Int, radius * os[1]) + 1
    ny = 2 * ceil(Int, radius * os[2]) + 1
    shape = (max(nx, 5), max(ny, 5))
    isodd(shape[1]) || (shape = (shape[1] + 1, shape[2]))
    isodd(shape[2]) || (shape = (shape[1], shape[2] + 1))
    origin = ((shape[1] + 1) / 2, (shape[2] + 1) / 2)

    # Build the initial ePSF directly from current star positions and fluxes.
    data = _stack_psf(
        image,
        stars,
        shape,
        origin,
        os;
        sigma_clip,
        sample_clip,
        smooth,
        recenter,
        clip_negative,
        badmask
    )
    psf = ImagePSF(data; origin, oversampling = os, fill_value)
    iterations = 0

    for iter in 1:maxiters
        # Refit each star against the current ePSF.
        iterations = iter
        oldx = [s.x for s in stars]
        oldy = [s.y for s in stars]
        for star in stars
            star.used || continue
            _fit_star!(
                star,
                psf,
                image;
                badmask,
                reweight,
                star_max_iter,
                fit_bkg,
                show_trace,
                kwargs...
            )
        end
        # Remove the arbitrary common centroid drift from the single-image fit.
        if anchor_centroids
            dxs = T[]
            dys = T[]
            for (k, star) in pairs(stars)
                star.used || continue
                push!(dxs, star.x - oldx[k])
                push!(dys, star.y - oldy[k])
            end
            if !isempty(dxs)
                dx_med = median(dxs)
                dy_med = median(dys)
                for star in stars
                    star.used || continue
                    star.x -= dx_med
                    star.y -= dy_med
                end
            end
        end
        count(s -> s.used, stars) ≥ 2 || throw(ArgumentError("too few stars survived empirical PSF fitting"))
        # Restack the ePSF using the updated per-star parameters.
        data = _stack_psf(
            image,
            stars,
            shape,
            origin,
            os;
            sigma_clip,
            sample_clip,
            smooth,
            recenter,
            clip_negative,
            badmask
        )
        psf = ImagePSF(data; origin, oversampling = os, fill_value)

        # Stop once the largest accepted centroid update is small.
        max_shift = zero(T)
        for (k, star) in pairs(stars)
            star.used || continue
            max_shift = max(max_shift, hypot(star.x - oldx[k], star.y - oldy[k]))
        end
        max_shift ≤ T(centroid_tol) && break
    end

    # Return the PSF plus per-star diagnostics for downstream filtering.
    result = ImagePSFBuildResult(
        psf,
        T[s.x for s in stars],
        T[s.y for s in stars],
        T[s.flux for s in stars],
        T[s.bkg for s in stars],
        BitVector(s.used for s in stars),
        BitVector(s.converged for s in stars),
        iterations,
        T[s.cost for s in stars],
    )
    return psf, result
end

"""
    fit(ImagePSF, image, x, y; fit_rad, kwargs...) -> (psf, result)

Build a single empirical ePSF from stars in `image`, starting from initial
detector-pixel coordinates `x` and `y`. Pixels whose full square area lies
inside `fit_rad` of the initial center are used for each stellar cutout.
"""
function fit(::Type{ImagePSF}, image::AbstractMatrix, x, y; fit_rad::Real, drop_edge::Bool = true, kwargs...)
    # Convert catalog positions into fit cutouts before entering the builder.
    stars = _stars_from_xy(image, x, y, fit_rad; drop_edge)
    return _fit_imagepsf(image, stars; psf_radius = fit_rad, kwargs...)
end

"""
    fit(ImagePSF, image, inds; x=nothing, y=nothing, kwargs...) -> (psf, result)

Build a single empirical ePSF from explicit cutout ranges. `inds` may be one
cutout `(xrange, yrange)` or an iterable of such cutouts. If `x` and `y` are
omitted, each initial center is the midpoint of its cutout.
"""
function fit(::Type{ImagePSF}, image::AbstractMatrix, inds; x = nothing, y = nothing, kwargs...)
    # Use caller-provided cutout ranges directly.
    stars = _stars_from_inds(image, inds; x, y)
    return _fit_imagepsf(image, stars; kwargs...)
end

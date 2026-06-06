# ==============================================================================
# Empirical ePSF Builder — Anderson & King (2000) iterative residual-stacking
#
# [`build_epsf`](@ref) implements the single-image ePSF construction method of
# [Anderson2000](@cite): extract star cutouts, project normalized pixels onto an
# oversampled grid, robustly combine, fill holes, smooth, recenter, normalize,
# and iteratively refit star parameters against the evolving ePSF until centroids
# converge.
# ==============================================================================

# ------------------------------------------------------------------------------
# Smoothing kernel — Anderson & King (2000) Eq. 8
# ------------------------------------------------------------------------------

"""The quartic smoothing kernel of [`Anderson2000`](@cite) (Eq. 8) used by default when `smooth=true` in `fit(ImagePSF, ...)`."""
const _QUARTIC_SMOOTHING_KERNEL = (
    (0.041632, -0.080816, 0.078368, -0.080816, 0.041632),
    (-0.080816, -0.019592, 0.200816, -0.019592, -0.080816),
    (0.078368, 0.200816, 0.441632, 0.200816, 0.078368),
    (-0.080816, -0.019592, 0.200816, -0.019592, -0.080816),
    (0.041632, -0.080816, 0.078368, -0.080816, 0.041632),
)

# ------------------------------------------------------------------------------
# Utility helpers
# ------------------------------------------------------------------------------

@inline function _round_half_away(x)
    return x ≥ zero(x) ? floor(Int, x + oftype(x, 0.5)) : ceil(Int, x - oftype(x, 0.5))
end

@inline function _pixel_wholly_inside(i, j, x, y, r)
    return (abs(i - x) + 0.5)^2 + (abs(j - y) + 0.5)^2 ≤ r^2
end

# ------------------------------------------------------------------------------
# EmpiricalStar — per-star cutout metadata
# ------------------------------------------------------------------------------

"""
    EmpiricalStar{T}

Internal per-star record for the ePSF builder. Holds the pixel cutout, current
centroid / flux / background estimates, and convergence diagnostics. Immutable;
updated via `ConstructionBase.setproperties`.
"""
struct EmpiricalStar{T}
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

ConstructionBase.getproperties(star::EmpiricalStar) =
    (x = star.x, y = star.y, flux = star.flux, bkg = star.bkg)

function ConstructionBase.setproperties(star::EmpiricalStar{T}, patch::NamedTuple) where {T}
    x = haskey(patch, :x) ? T(patch.x) : star.x
    y = haskey(patch, :y) ? T(patch.y) : star.y
    flux = haskey(patch, :flux) ? T(patch.flux) : star.flux
    bkg = haskey(patch, :bkg) ? T(patch.bkg) : star.bkg
    used = haskey(patch, :used) ? patch.used : star.used
    converged = haskey(patch, :converged) ? patch.converged : star.converged
    cost = haskey(patch, :cost) ? T(patch.cost) : star.cost
    return EmpiricalStar{T}(star.inds, star.pixels, x, y, flux, bkg, used, converged, cost)
end

# ------------------------------------------------------------------------------
# BuilderState — grid geometry and algorithm configuration
# ------------------------------------------------------------------------------

"""
    BuilderState{T}

Immutable configuration bundle carrying the oversampled ePSF grid geometry and
all algorithm toggles for one `build_epsf` run.
"""
struct BuilderState{T}
    shape::Tuple{Int, Int}
    origin::Tuple{T, T}
    oversampling::Tuple{Int, Int}
    sigma_clip::T
    sample_clip::T
    smooth::Bool
    recenter::Bool
    clip_negative::Bool
    anchor_centroids::Bool
    fill_value::T
end

function BuilderState(
        oversampling, psf_radius, stars;
        sigma_clip, sample_clip, smooth, recenter, clip_negative, anchor_centroids, fill_value
    )
    T = typeof(stars[1].x)
    os = _as_oversampling(oversampling)
    # Choose an odd supersampled PSF grid with a central origin.
    radius = if isnothing(psf_radius)
        maximum(stars) do star
            maximum(max(abs(idx[1] - star.x), abs(idx[2] - star.y)) for idx in star.pixels)
        end
    else
        T(psf_radius)
    end
    nx = 2 * ceil(Int, radius * os[1]) + 1
    ny = 2 * ceil(Int, radius * os[2]) + 1
    shape = (max(nx, 5), max(ny, 5))
    isodd(shape[1]) || (shape = (shape[1] + 1, shape[2]))
    isodd(shape[2]) || (shape = (shape[1], shape[2] + 1))
    # Central pixel is the geometric origin so (shape+1)/2 is exact.
    origin = (T((shape[1] + 1) / 2), T((shape[2] + 1) / 2))
    return BuilderState{T}(
        shape, origin, os,
        T(sigma_clip), T(sample_clip), smooth, recenter,
        clip_negative, anchor_centroids, T(fill_value)
    )
end

# ==============================================================================
# Phase: Star extraction
# ==============================================================================

"""
    estimate_local_sky(image, inds)

Use finite edge pixels of a cutout as a simple local sky estimate.
"""
function estimate_local_sky(image, inds)
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

"""
    estimate_initial_star_params(star, image) -> EmpiricalStar

Estimate a star's background from the cutout edge and its flux from a
winsorized (98th-percentile-capped) aperture sum of background-subtracted
pixels. Returns a new `EmpiricalStar` with updated `bkg` and `flux`.
"""
function estimate_initial_star_params(star::EmpiricalStar{T}, image) where {T}
    bkg = estimate_local_sky(image, star.inds)
    vals = T[]
    for idx in star.pixels
        v = image[idx]
        isfinite(v) && push!(vals, max(zero(T), T(v - bkg)))
    end
    # Winsorize the brightest pixels so hot defects do not set the aperture flux.
    sort!(vals)
    cap = isempty(vals) ? zero(T) : vals[clamp(ceil(Int, 0.98 * length(vals)), 1, length(vals))]
    flux = sum(v -> min(v, cap), vals; init = zero(T))
    # Fall back to a peak estimate if the aperture sum is unusable.
    if !(isfinite(flux) && flux > eps(T))
        flux = max(eps(T), T(maximum(image[star.inds...] .- bkg)))
    end
    return setproperties(star, (bkg = T(bkg), flux = T(flux)))
end

function _normalize_cutout_inds(inds)
    # Treat one `(xrange, yrange)` tuple as a single cutout, not two cutouts.
    if inds isa Tuple && length(inds) == 2 && inds[1] isa AbstractUnitRange && inds[2] isa AbstractUnitRange
        return (inds,)
    end
    return Tuple(inds)
end

"""
    extract_stars(image, x, y, fit_rad; drop_edge=false) -> Vector{EmpiricalStar}

Build per-star cutouts from an image and initial detector-pixel coordinates.
Pixels whose full square area lies inside `fit_rad` of the initial center are
used for each star. Each cutout gets initial centroid, flux, and background
estimates.
"""
function extract_stars(image, x, y, fit_rad; drop_edge::Bool)
    xs = collect(x)
    ys = collect(y)
    length(xs) == length(ys) || throw(ArgumentError("`x` and `y` must have the same length"))
    T = promote_type(eltype(image), eltype(xs), eltype(ys), typeof(fit_rad))
    T = T <: Integer ? Float64 : float(T)
    ax, ay = axes(image)
    stars = EmpiricalStar{T}[]
    for k in eachindex(xs, ys)
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
        star = EmpiricalStar((xlo:xhi, ylo:yhi), pix, xk, yk, one(T), zero(T), true, false, T(Inf))
        star = estimate_initial_star_params(star, image)
        push!(stars, star)
    end
    isempty(stars) && throw(ArgumentError("no usable stars remain after cutout extraction"))
    return stars
end

"""
    extract_stars(image, inds; x=nothing, y=nothing) -> Vector{EmpiricalStar}

Build per-star cutouts from explicit pixel ranges. `inds` may be one cutout
`(xrange, yrange)` or an iterable of such cutouts. If `x` and `y` are omitted,
each initial center is the midpoint of its cutout.
"""
function extract_stars(image, inds; x = nothing, y = nothing)
    cutouts = _normalize_cutout_inds(inds)
    n = length(cutouts)
    if !isnothing(x) || !isnothing(y)
        (isnothing(x) || isnothing(y)) && throw(ArgumentError("`x` and `y` must be provided together"))
        length(x) == n && length(y) == n || throw(ArgumentError("`x` and `y` must match `inds` length"))
    end
    T = promote_type(eltype(image), isnothing(x) ? Float64 : eltype(x), isnothing(y) ? Float64 : eltype(y))
    T = T <: Integer ? Float64 : float(T)
    stars = EmpiricalStar{T}[]
    for k in 1:n
        c = cutouts[k]
        length(c) == 2 || throw(ArgumentError("each cutout must be a tuple of two ranges"))
        xr = Int(first(c[1])):Int(last(c[1]))
        yr = Int(first(c[2])):Int(last(c[2]))
        first(axes(image, 1)) ≤ first(xr) && last(xr) ≤ last(axes(image, 1)) ||
            throw(ArgumentError("x range for cutout $k is outside the image"))
        first(axes(image, 2)) ≤ first(yr) && last(yr) ≤ last(axes(image, 2)) ||
            throw(ArgumentError("y range for cutout $k is outside the image"))
        xk = isnothing(x) ? T((first(xr) + last(xr)) / 2) : T(x[k])
        yk = isnothing(y) ? T((first(yr) + last(yr)) / 2) : T(y[k])
        pix = vec(collect(CartesianIndices((xr, yr))))
        star = EmpiricalStar((xr, yr), pix, xk, yk, one(T), zero(T), true, false, T(Inf))
        star = estimate_initial_star_params(star, image)
        push!(stars, star)
    end
    isempty(stars) && throw(ArgumentError("`inds` must contain at least one cutout"))
    return stars
end

# ==============================================================================
# Phase: Grid infill utilities
# ==============================================================================

function _robust_location(vals::Vector{T}, sigma_clip) where {T}
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

"""
    nearest_finite_samples(data, i, j; along_x::Bool)

Search outward from `(i, j)` along one axis until four finite samples are found.
"""
function nearest_finite_samples(data, i, j; along_x::Bool)
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

"""
    has_finite_support_on_both_sides(data, i, j; along_x::Bool)

Check whether `data` has finite samples on both sides of `(i, j)` along one
axis, so infill can interpolate rather than extrapolate.
"""
function has_finite_support_on_both_sides(data, i, j; along_x::Bool)
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

"""
    cubic_lagrange_interpolate(xs, ys, x0)

Evaluate the cubic Lagrange interpolant through four arbitrary `(xs, ys)`
samples at `x0`. Used during grid hole-filling when the finite support does not
form a contiguous stencil for bicubic interpolation.
"""
function cubic_lagrange_interpolate(xs, ys, x0)
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
    fill_grid_holes!(data; maxiter::Int = 6))

Fill non-finite holes in an empirical ePSF stack in place. Interior holes are
filled first with bracketed separable cubic interpolation, unsupported border
holes then use a local exponential radial extrapolation
(see 4.2.1 of [Anderson2000](@citet)) and any remaining holes fall back to a
local median if possible, and zero fill otherwise. `stack_epsf_grid` uses
this after robustly combining stellar samples into the *initial* oversampled grid
and before smoothing, recentering, and normalization. `maxiter` controls the
maximum number of passes for the bicubic and median stages; six is a reasonable
default to allow infill growth without risking infinite loops.
"""
function fill_grid_holes!(data; maxiter::Int = 6)
    T = eltype(data)
    nx, ny = size(data)
    all(isfinite, data) && return data

    fillable = falses(nx, ny)
    for j in 1:ny, i in 1:nx
        if !isfinite(data[i, j])
            fillable[i, j] = has_finite_support_on_both_sides(data, i, j; along_x = true) &&
                has_finite_support_on_both_sides(data, i, j; along_x = false)
        end
    end
    if count(fillable) / (nx * ny) > 0.1
        @warn "more than 10% of the ePSF grid is missing; infill may be unreliable"
    end
    # First: fill interior holes with separable cubic interpolation where possible
    for _ in 1:maxiter
        changed = false
        old = copy(data)
        for j in 1:ny, i in 1:nx
            isfinite(old[i, j]) && continue
            fillable[i, j] || continue
            yinds = nearest_finite_samples(old, i, j; along_x = false)
            length(yinds) == 4 || continue
            any(yy -> yy < j, yinds) && any(yy -> yy > j, yinds) || continue
            row_vals = T[]
            good_y = Int[]
            for yy in yinds
                has_finite_support_on_both_sides(old, i, yy; along_x = true) || continue
                xinds = nearest_finite_samples(old, i, yy; along_x = true)
                length(xinds) == 4 || continue
                any(xx -> xx < i, xinds) && any(xx -> xx > i, xinds) || continue
                vals = T[old[xx, yy] for xx in xinds]
                push!(row_vals, cubic_lagrange_interpolate(xinds, vals, i))
                push!(good_y, yy)
            end
            if length(row_vals) == 4
                data[i, j] = cubic_lagrange_interpolate(good_y, row_vals, j)
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
            v = old[ii, jj]
            if isfinite(v) && v > zero(T)
                push!(rs, hypot(T(ii) - cx, T(jj) - cy))
                push!(logs, log(T(v)))
                push!(vals, T(v))
            end
        end
        length(rs) ≥ 2 || continue

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
    for idx in eachindex(data)
        isfinite(data[idx]) || (data[idx] = zero(T))
    end
    return data
end

# ==============================================================================
# Phase: ePSF grid stacking pipeline
# ==============================================================================

"""
    project_star_pixels_to_grid(image, stars, state, badmask) -> Vector{Vector{T}}

Project each valid, normalized stellar pixel onto the oversampled ePSF grid,
returning per-cell sample lists.
"""
function project_star_pixels_to_grid(image, stars, state::BuilderState{T}, badmask) where {T}
    nx, ny = state.shape
    cells = [T[] for _ in 1:(nx * ny)]
    sx, sy = state.oversampling
    ox, oy = state.origin
    for star in stars
        star.used || continue
        isfinite(star.flux) && star.flux > eps(T) || continue
        for idx in star.pixels
            isnothing(badmask) || !badmask[idx] || continue
            v = image[idx]
            isfinite(v) || continue
            gx = _round_half_away(T(ox) + T(sx) * (T(idx[1]) - star.x))
            gy = _round_half_away(T(oy) + T(sy) * (T(idx[2]) - star.y))
            1 ≤ gx ≤ nx && 1 ≤ gy ≤ ny || continue
            sample = T(v - star.bkg) / star.flux
            if abs(sample) > state.sample_clip
                continue
            end
            push!(cells[gx + (gy - 1) * nx], sample)
        end
    end
    return cells
end

"""
    robust_combine_grid_cells(cells, state) -> AbstractMatrix

Apply sigma-clipped median combination to each ePSF cell, returning a grid
with NaN markers for unfilled cells.
"""
function robust_combine_grid_cells(cells, state::BuilderState{T}) where {T}
    nx, ny = state.shape
    data = fill(T(NaN), nx, ny)
    for gy in 1:ny, gx in 1:nx
        vals = cells[gx + (gy - 1) * nx]
        isempty(vals) || (data[gx, gy] = _robust_location(vals, state.sigma_clip))
    end
    return data
end

function smooth_grid_quartic!(data::AbstractMatrix{T}) where {T}
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

function recenter_grid_to_origin!(data, origin, oversampling)
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
    out = similar(data)
    for j in axes(data, 2), i in axes(data, 1)
        out[i, j], _, _ = bicubic_interpolate(data, T(i) + dx, T(j) + dy; fill_value = zero(T))
    end
    return out
end

function normalize_grid_to_oversampling_area!(data, oversampling; clip_negative::Bool)
    T = eltype(data)
    if clip_negative
        for idx in eachindex(data)
            data[idx] = max(zero(T), data[idx])
        end
    end
    s = sum(data)
    if !(isfinite(s) && s > eps(T))
        throw(ArgumentError("empirical PSF stack has non-positive finite sum"))
    end
    data .*= T(oversampling[1] * oversampling[2]) / s
    return data
end

"""
    stack_epsf_grid(image, stars, state; badmask=nothing)

Run the full ePSF stacking pipeline: project star pixels → robustly combine →
fill holes → optionally smooth → optionally recenter → normalize to the
oversampling area convention.
"""
function stack_epsf_grid(image, stars, state::BuilderState; badmask = nothing)
    cells = project_star_pixels_to_grid(image, stars, state, badmask)
    data = robust_combine_grid_cells(cells, state)
    fill_grid_holes!(data)
    if state.smooth === true || state.smooth === :quartic
        data = smooth_grid_quartic!(data)
    elseif !(state.smooth === false || isnothing(state.smooth))
        throw(ArgumentError("`smooth` must be true, false, `:quartic`, or `nothing`"))
    end
    state.recenter && (data = recenter_grid_to_origin!(data, state.origin, state.oversampling))
    normalize_grid_to_oversampling_area!(data, state.oversampling; clip_negative = state.clip_negative)
    return data
end

# ==============================================================================
# Phase: Per-star fitting against the current ePSF
# ==============================================================================

function finite_unmasked_pixels(star, image, badmask)
    pix = CartesianIndex{2}[]
    for idx in star.pixels
        isnothing(badmask) || !badmask[idx] || continue
        isfinite(image[idx]) && push!(pix, idx)
    end
    return pix
end

"""
    fit_star_against_epsf(star, psf, image; kwargs...) -> EmpiricalStar

Fit a single star's centroid, flux, and optionally background against the
current ePSF using LM/IRLS. Returns a new `EmpiricalStar` with updated
parameters. If the fit fails, returns the star with `used=false`.
"""
function fit_star_against_epsf(
        star::EmpiricalStar{T},
        psf::ImagePSF,
        image;
        badmask = nothing,
        reweight = TukeyLoss(),
        star_max_iter::Integer = 100,
        fit_bkg::Bool = false,
        show_trace::Bool = false,
        kwargs...
    ) where {T}

    pix = finite_unmasked_pixels(star, image, badmask)
    nparams = fit_bkg ? 4 : 3
    if length(pix) <= nparams
        return setproperties(star, (used = false, converged = false))
    end
    x0 = fit_bkg ? T[star.x, star.y, star.flux, star.bkg] : T[star.x, star.y, star.flux]
    star_bkg = star.bkg  # capture before closure for fixed-bkg mode

    function accum!(A::AbstractMatrix{FT}, b::AbstractVector{FT}, residuals::AbstractVector{FT}, x::AbstractVector{FT}, weights) where {FT}
        fill!(A, zero(FT))
        fill!(b, zero(FT))
        cost = zero(FT)
        model = ImagePSF(
            psf.data; x = x[1], y = x[2], flux = x[3],
            bkg = fit_bkg ? x[4] : star_bkg, origin = psf.origin,
            oversampling = psf.oversampling, fill_value = psf.fill_value
        )
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

    problem = LMProblem(x0, length(pix), accum!)
    result = lm_irls(
        problem;
        max_iter = star_max_iter,
        reweight,
        show_trace,
        kwargs...
    )
    p = result.minimizer
    if all(isfinite, p) && p[3] > eps(T)
        return setproperties(
            star, (
                x = T(p[1]), y = T(p[2]), flux = T(p[3]),
                bkg = fit_bkg ? T(p[4]) : star.bkg,
                converged = result.converged,
                cost = T(result.minimum),
                used = result.converged,
            )
        )
    else
        return setproperties(star, (used = false, converged = false, cost = T(Inf)))
    end
end

# ==============================================================================
# Phase: Iterative ePSF builder
# ==============================================================================

function fit_all_stars(
        stars, psf, image; badmask = nothing, reweight = TukeyLoss(),
        star_max_iter::Integer = 100, fit_bkg::Bool = false, show_trace::Bool = false, kwargs...
    )
    for k in eachindex(stars)
        stars[k].used || continue
        stars[k] = fit_star_against_epsf(
            stars[k], psf, image;
            badmask, reweight, star_max_iter, fit_bkg, show_trace, kwargs...
        )
    end
    return stars
end

function remove_centroid_drift(stars, old_centroids)
    T = typeof(stars[1].x)
    dxs = T[]
    dys = T[]
    for k in eachindex(stars)
        stars[k].used || continue
        push!(dxs, stars[k].x - old_centroids[k][1])
        push!(dys, stars[k].y - old_centroids[k][2])
    end
    if !isempty(dxs)
        dx_med = median(dxs)
        dy_med = median(dys)
        for k in eachindex(stars)
            stars[k].used || continue
            stars[k] = setproperties(stars[k], (x = stars[k].x - dx_med, y = stars[k].y - dy_med))
        end
    end
    return stars
end

function build_result(stars, psf, iterations)
    T = typeof(stars[1].x)
    return ImagePSFBuildResult(
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
end

"""
    build_epsf(image, stars; kwargs...) -> (psf, result)

Build a single empirical ePSF from per-star cutouts using the Anderson & King
(2000) iterative residual-stacking method ([Anderson2000](@cite)).

The algorithm:
1. Projects normalized star pixels onto an oversampled grid.
2. Robustly combines samples (sigma-clipped median), fills holes, optionally
   smooths and recenters, then normalizes.
3. Fits each star's centroid and flux against the current ePSF using LM/IRLS.
4. Removes the median centroid drift to break the PSF/centroid degeneracy.
5. Restacks the ePSF grid from the updated star parameters.
6. Repeats 3–5 until centroid shifts fall below `centroid_tol`.
"""
function build_epsf(
        image::AbstractMatrix,
        stars::Vector{EmpiricalStar{T}};
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
    if !isnothing(badmask)
        size(badmask) == size(image) ||
            throw(ArgumentError("`badmask` must have the same size as `image`"))
    end
    state = BuilderState(
        oversampling, psf_radius, stars;
        sigma_clip, sample_clip, smooth, recenter,
        clip_negative, anchor_centroids, fill_value
    )

    # Phase 1: Initial ePSF stack from current star estimates.
    data = stack_epsf_grid(image, stars, state; badmask)
    psf = ImagePSF(
        data; origin = state.origin, oversampling = state.oversampling,
        fill_value = state.fill_value
    )

    # Phase 2–N: Iteratively refit stars and restack.
    iterations = 0
    for iter in 1:maxiters
        iterations = iter
        old_centroids = [(s.x, s.y) for s in stars]
        stars = fit_all_stars(
            stars, psf, image;
            badmask, reweight, star_max_iter, fit_bkg, show_trace, kwargs...
        )
        # Capture fitted centroids before anchoring so the convergence
        # check reflects the actual fit improvement, not the anchoring step.
        fitted_centroids = [(s.x, s.y) for s in stars]
        if anchor_centroids
            stars = remove_centroid_drift(stars, old_centroids)
        end
        count(s -> s.used, stars) ≥ 2 ||
            throw(ArgumentError("too few stars survived empirical PSF fitting"))
        data = stack_epsf_grid(image, stars, state; badmask)
        psf = ImagePSF(
            data; origin = state.origin, oversampling = state.oversampling,
            fill_value = state.fill_value
        )
        # Convergence uses pre-anchoring shifts so the loop does not
        # exit before the ePSF has a chance to improve.
        max_shift = zero(T)
        for k in eachindex(stars)
            stars[k].used || continue
            max_shift = max(max_shift, hypot(
                fitted_centroids[k][1] - old_centroids[k][1],
                fitted_centroids[k][2] - old_centroids[k][2],
            ))
        end
        max_shift ≤ T(centroid_tol) && break
    end
    return psf, build_result(stars, psf, iterations)
end

# ==============================================================================
# Public API
# ==============================================================================

"""
    fit(ImagePSF, image, x, y; fit_rad, kwargs...) -> (psf, result)

Build a single empirical ePSF from stars in `image`, starting from initial
detector-pixel coordinates `x` and `y`. Pixels whose full square area lies
inside `fit_rad` of the initial center are used for each stellar cutout.
"""
function fit(
        ::Type{ImagePSF}, image::AbstractMatrix, x, y;
        fit_rad::Real, drop_edge::Bool = true, kwargs...
    )
    stars = extract_stars(image, x, y, fit_rad; drop_edge)
    return build_epsf(image, stars; psf_radius = fit_rad, kwargs...)
end

"""
    fit(ImagePSF, image, inds; x=nothing, y=nothing, kwargs...) -> (psf, result)

Build a single empirical ePSF from explicit cutout ranges. `inds` may be one
cutout `(xrange, yrange)` or an iterable of such cutouts. If `x` and `y` are
omitted, each initial center is the midpoint of its cutout.
"""
function fit(
        ::Type{ImagePSF}, image::AbstractMatrix, inds;
        x = nothing, y = nothing, kwargs...
    )
    stars = extract_stars(image, inds; x, y)
    return build_epsf(image, stars; kwargs...)
end

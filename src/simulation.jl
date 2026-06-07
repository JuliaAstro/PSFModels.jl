function _sample_flux(rng::Random.AbstractRNG, flux, distribution, power)
    # A scalar flux means every simulated source has the same brightness.
    if flux isa Number
        return float(flux)
    elseif flux isa Tuple || flux isa AbstractVector
        # A two-value flux range is sampled from the requested distribution.
        length(flux) == 2 || throw(ArgumentError("flux range must have two values"))
        lo, hi = float(flux[1]), float(flux[2])
        lo > 0 && hi ≥ lo || throw(ArgumentError("flux range must be positive and ordered"))
        u = rand(rng)
        if distribution === :uniform
            return lo + u * (hi - lo)
        elseif distribution === :loguniform
            return exp(log(lo) + u * (log(hi) - log(lo)))
        elseif distribution === :powerlaw
            α = float(power)
            if abs(α - 1) < eps(typeof(α)) * 10
                return exp(log(lo) + u * (log(hi) - log(lo)))
            end
            q = 1 - α
            return (lo^q + u * (hi^q - lo^q))^(1 / q)
        else
            throw(ArgumentError("unsupported flux distribution `$distribution`"))
        end
    else
        throw(ArgumentError("`flux` must be a scalar or two-value range"))
    end
end

@doc raw"""
    flux_for_snr(model; snr, background, read_noise=0, gain=1)

Approximate source flux, in the same data units as `background` and
`read_noise` (for example, ADU), required for matched-filter `snr` for a PSF
model on a flat background. `gain` is in electrons per data unit. The estimate
uses [`effective_area(model)`](@ref PSFModels.effective_area), includes source
shot noise, and a background/read-noise variance per pixel of
`background / gain + read_noise^2`. The returned flux is in data units; the
corresponding source flux in electrons is `F * gain`:

```math
\\mathrm{snr} \\approx \\frac{F}{\\sqrt{F / g + A_\\mathrm{eff}\\,\\sigma_b^2}},
\\quad
\\sigma_b^2 = \\max\\left(0, \\frac{\\mathrm{background}}{g} + \\mathrm{read\\_noise}^2\\right)
```

Solving the quadratic for the positive flux gives

```math
F \approx \frac{1}{2}\left[
\frac{\mathrm{snr}^2}{g} +
\sqrt{\left(\frac{\mathrm{snr}^2}{g}\right)^2
+ 4\,\mathrm{snr}^2 A_\mathrm{eff}\sigma_b^2}
\right],
\quad
A_\mathrm{eff} = \mathrm{effective\_area}(\mathrm{model})
```
"""
function flux_for_snr(
        model::AbstractPSFModel;
        snr::Real,
        background::Real,
        read_noise::Real = 0,
        gain::Real = 1
    )
    # Validate inputs before solving the SNR-flux quadratic.
    snr ≥ 0 || throw(ArgumentError("`snr` must be non-negative"))
    gain > 0 || throw(ArgumentError("`gain` must be positive"))
    snr, background, read_noise, gain = float.((snr, background, read_noise, gain))

    # Combine the background/read-noise variance and effective PSF area.
    σb2 = max(zero(background), background) / gain + read_noise^2
    background_variance = effective_area(model) * σb2

    # Include source shot noise by taking the positive root of the quadratic.
    source_noise_term = snr^2 / gain
    return (source_noise_term + sqrt(source_noise_term^2 + 4 * snr^2 * background_variance)) / 2
end

function _flux_from_snr_spec(rng, model, snr, background, read_noise, gain)
    # Collapse matrix backgrounds to a representative scalar for source draws.
    bg = background isa AbstractMatrix ? median(background) : background
    if snr isa Number
        return flux_for_snr(model; snr, background = bg, read_noise, gain)
    elseif snr isa Tuple || snr isa AbstractVector
        # Allow source-to-source SNR variation over a uniform range.
        length(snr) == 2 || throw(ArgumentError("`snr` range must have two values"))
        sampled = snr[1] + rand(rng) * (snr[2] - snr[1])
        return flux_for_snr(model; snr = sampled, background = bg, read_noise, gain)
    else
        throw(ArgumentError("`snr` must be a scalar or two-value range"))
    end
end

"""
    simulate_sources(shape, n_sources; kwargs...) -> NamedTuple

Generate random source positions and fluxes for an image with dimensions
`shape == (nx, ny)`.

Source fluxes can be specified in one of two mutually exclusive ways:

**Branch 1 — Flux-based (default):** `flux`, `flux_distribution`, `flux_power`
are used to draw fluxes directly; `psf`, `background`, `read_noise`, and
`gain` are ignored.

- `flux`: scalar or `(lo, hi)` range. Defaults to `(100, 1000)`.
- `flux_distribution`: `:uniform`, `:loguniform`, or `:powerlaw`.
- `flux_power`: power-law exponent for `:powerlaw` (`dN/dF ∝ F^-flux_power`).

**Branch 2 — SNR-based (when `snr` is set):** fluxes are derived from a
requested signal-to-noise ratio via [`flux_for_snr`](@ref), which requires
the PSF model and noise properties. `flux` and `flux_distribution` are
ignored in this branch.

- `snr`: scalar SNR or `(lo, hi)` range. If a scalar, every source is drawn
  at the flux required for that SNR. If a two-value range, each source gets a
  uniform random SNR in `[lo, hi]` before converting to flux.
- `psf`: an [`AbstractPSFModel`](@ref). Required when `snr` is provided in
  order to compute the effective area for flux conversion.
- `background`: scalar background level or a 2-D `AbstractMatrix`
  (image-sized). If a matrix, its median is used as a representative scalar.
  Defaults to `0`.
- `read_noise`: Gaussian read-noise standard deviation in data units (e.g.,
  ADU). Defaults to `0`.
- `gain`: CCD gain in electrons per data unit. Must be positive. Defaults to `1`.

**Common keywords (used in both branches):**

- `min_separation`: minimum source-center separation in pixels.
- `border`: scalar or `(xborder, yborder)` excluded image border.
"""
function simulate_sources(
        shape::Tuple{<:Integer, <:Integer},
        n_sources::Integer;
        flux = (100.0, 1000.0),
        flux_distribution = :uniform,
        flux_power::Real = 2,
        min_separation::Real = 0,
        border = 0,
        rng::Random.AbstractRNG = Random.default_rng(),
        snr = nothing,
        psf::Union{Nothing, AbstractPSFModel} = nothing,
        background::Union{Real, AbstractMatrix} = 0,
        read_noise::Real = 0,
        gain::Real = 1,
        max_attempts::Integer = max(10_000, 200 * n_sources)
    )
    # Define the valid placement region after excluding the image border.
    n_sources ≥ 0 || throw(ArgumentError("`n_sources` must be non-negative"))
    # Treat scalar inputs as identical values on both axes.
    bx, by = if length(border) == 1
        (border, border)
    elseif length(border) == 2
        (border[1], border[2])
    else
        throw(ArgumentError("`border` must be a scalar or length-2 tuple/vector"))
    end
    xmin, xmax = 1 + float(bx), float(shape[1]) - float(bx)
    ymin, ymax = 1 + float(by), float(shape[2]) - float(by)
    xmin ≤ xmax && ymin ≤ ymax || throw(ArgumentError("border leaves no valid source-placement area"))
    minsep2 = float(min_separation)^2
    xs = Float64[]
    ys = Float64[]
    fs = Float64[]
    attempts = 0
    while length(xs) < n_sources && attempts < max_attempts
        # Draw a candidate position and reject it if it is too close to earlier sources.
        attempts += 1
        x = xmin + rand(rng) * (xmax - xmin)
        y = ymin + rand(rng) * (ymax - ymin)
        separated = true
        for k in eachindex(xs)
            if (x - xs[k])^2 + (y - ys[k])^2 < minsep2
                separated = false
                break
            end
        end
        separated || continue
        # Draw either flux directly or infer it from a requested SNR.
        f = if isnothing(snr)
            _sample_flux(rng, flux, flux_distribution, flux_power)
        else
            isnothing(psf) && throw(ArgumentError("`psf` is required when generating fluxes from `snr`"))
            _flux_from_snr_spec(rng, psf, snr, background, read_noise, gain)
        end
        push!(xs, x)
        push!(ys, y)
        push!(fs, f)
    end
    if attempts >= max_attempts
        @warn "Reached maximum attempts ($max_attempts) with only $(length(xs)) sources generated; consider reducing `min_separation` or increasing `max_attempts`. Returning the $(length(xs)) sources that were generated."
    end
    return (id = collect(1:length(xs)), x = xs, y = ys, flux = fs)
end

function _source_vectors(sources)
    # NamedTuple/table-like inputs expose x, y, and flux directly.
    if hasproperty(sources, :x) && hasproperty(sources, :y) && hasproperty(sources, :flux)
        return sources.x, sources.y, sources.flux
    end
    # Fall back to iterating row-like objects with x/y/flux properties.
    xs = Float64[]
    ys = Float64[]
    fs = Float64[]
    for src in sources
        push!(xs, getproperty(src, :x))
        push!(ys, getproperty(src, :y))
        push!(fs, getproperty(src, :flux))
    end
    return xs, ys, fs
end

function _source_ranges(model, x, y, model_radius, image)
    # Use the model extent unless the caller supplies a fixed render radius.
    if isnothing(model_radius)
        (xlo, xhi), (ylo, yhi) = extent(model)
    else
        xlo, xhi = x - model_radius, x + model_radius
        ylo, yhi = y - model_radius, y + model_radius
    end
    # Clip the render footprint to the output image.
    xr = max(first(axes(image, 1)), floor(Int, xlo)):min(last(axes(image, 1)), ceil(Int, xhi))
    yr = max(first(axes(image, 2)), floor(Int, ylo)):min(last(axes(image, 2)), ceil(Int, yhi))
    return xr, yr
end

"""
    render_sources!(image, model, sources; model_radius=nothing)

Add sources to `image` in place using `model`. `sources` must provide `x`, `y`,
and `flux` fields or columns. The model background is set to zero while each
source is rendered.
"""
function render_sources!(
        image::AbstractMatrix,
        model::AbstractPSFModel,
        sources;
        model_radius = nothing
    )
    # Normalize source containers to coordinate and flux vectors.
    xs, ys, fs = _source_vectors(sources)
    length(xs) == length(ys) == length(fs) || throw(ArgumentError("source x/y/flux lengths must match"))
    for k in eachindex(xs, ys, fs)
        # Render each source with zero model background; image background is added separately.
        m = setproperties(model, (x = xs[k], y = ys[k], flux = fs[k], bkg = zero(float(eltype(image)))))
        xr, yr = _source_ranges(m, xs[k], ys[k], model_radius, image)
        for i in xr, j in yr
            image[i, j] += evaluate(m, i, j)
        end
    end
    return image
end

function _rand_poisson(rng::Random.AbstractRNG, λ::Real)
    # Use exact inversion for small means.
    λf = float(max(λ, zero(λ)))
    if λf < 30
        L = exp(-λf)
        k = 0
        p = 1.0
        while p > L
            k += 1
            p *= rand(rng)
        end
        return k - 1
    else
        # Use a Gaussian approximation for large means to avoid slow loops.
        return max(0, round(Int, λf + sqrt(λf) * randn(rng)))
    end
end

function rand_poisson(rng::Random.AbstractRNG, λ::Real)
    # Use exact inversion for small means and Gaussian approximation for large means.
    function _rand_poisson_knuth(rng::Random.AbstractRNG, λ::T) where T
        L = exp(-λ)
        k = 0
        p = one(T)

        while p > L
            k += 1
            p *= rand(rng, T)
        end

        return k - 1
    end

    function _rand_poisson_normal(rng::Random.AbstractRNG, λ::T) where T
        σ = sqrt(λ)

        while true
            x = round(Int, λ + σ * randn(rng, T))
            x >= 0 && return x
        end
    end

    λ < 0 && throw(DomainError(λ, "Poisson rate λ must be nonnegative"))
    λ == 0 && return 0

    if λ < 30
        return _rand_poisson_knuth(rng, float(λ))
    else
        return _rand_poisson_normal(rng, float(λ))
    end
end
rand_poisson(λ::Real) = rand_poisson(Random.default_rng(), λ)


"""
    add_noise!(image; noise=:poisson_gaussian, read_noise=0, gain=1, rng=Random.default_rng())

Apply in-place noise to an expected-count image. `:poisson` samples shot noise,
`:gaussian` adds read noise, and `:poisson_gaussian` applies both.
"""
function add_noise!(
        image::AbstractMatrix;
        noise = :poisson_gaussian,
        read_noise::Real = 0,
        gain::Real = 1,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    # Dispatch the requested noise model explicitly.
    gain > 0 || throw(ArgumentError("`gain` must be positive"))
    if isnothing(noise) || noise === :none
        return image
    elseif noise === :poisson
        for idx in eachindex(image)
            image[idx] = _rand_poisson(rng, max(zero(float(image[idx])), float(image[idx]) * gain)) / gain
        end
    elseif noise === :gaussian
        for idx in eachindex(image)
            image[idx] += read_noise * randn(rng)
        end
    elseif noise === :poisson_gaussian
        add_noise!(image; noise = :poisson, gain, rng)
        add_noise!(image; noise = :gaussian, read_noise, rng)
    else
        throw(ArgumentError("unsupported noise model `$noise`"))
    end
    return image
end

"""
    simulate_image(shape, psf, sources; background=0, noise=:poisson_gaussian, kwargs...)

Create an artificial image by rendering `sources` through `psf`, adding scalar or
matrix `background`, and optionally injecting Poisson/read noise.
"""
function simulate_image(
        shape::Tuple{<:Integer, <:Integer},
        psf::AbstractPSFModel,
        sources;
        background::Union{Real, AbstractMatrix} = 0,
        noise = :poisson_gaussian,
        read_noise::Real = 0,
        gain::Real = 1,
        model_radius = nothing,
        rng::Random.AbstractRNG = Random.default_rng()
    )
    # Allocate the image from either scalar or matrix background.
    T = promote_type(float(eltype(getproperties(psf))), background isa AbstractMatrix ? eltype(background) : typeof(background))
    T = T <: Integer ? Float64 : float(T)
    image = if background isa AbstractMatrix
        size(background) == shape || throw(ArgumentError("background matrix must match `shape`"))
        Matrix{T}(background)
    else
        fill(T(background), shape)
    end
    # Add source model values and then apply the requested noise model.
    render_sources!(image, psf, sources; model_radius)
    add_noise!(image; noise, read_noise, gain, rng)
    return image
end

"""
    simulate_image(shape, psf, n_sources; kwargs...) -> (image, sources)

Generate random sources with `simulate_sources`, render them through `psf`, and
return both the image and source table.
"""
function simulate_image(
        shape::Tuple{<:Integer, <:Integer},
        psf::AbstractPSFModel,
        n_sources::Integer;
        background::Union{Real, AbstractMatrix} = 0,
        noise = :poisson_gaussian,
        read_noise::Real = 0,
        gain::Real = 1,
        model_radius = nothing,
        rng::Random.AbstractRNG = Random.default_rng(),
        source_kwargs...
    )
    # Generate a source catalog with the same noise/background assumptions.
    sources = simulate_sources(
        shape, n_sources;
        rng, psf, background, read_noise, gain, source_kwargs...
    )
    # Render the generated catalog into an image.
    image = simulate_image(
        shape, psf, sources;
        background, noise, read_noise, gain, model_radius, rng
    )
    return image, sources
end

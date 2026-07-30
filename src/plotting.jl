# Plotting interface. The implementations live in package extensions:
#  - ext/PSFModelsPlotsExt.jl  (loaded with Plots.jl, via a RecipesBase recipe)
#  - ext/PSFModelsMakieExt.jl  (loaded with Makie.jl, via a Makie recipe)

export psfplot, psfplot!, psfplotview

"""
    psfplot(model, inds; kwargs...)
    psfplot(model, xs, ys; kwargs...)
    psfplot!([plt_or_ax,] model, args...; kwargs...)

Plot the PSF `model` as a heatmap evaluated over the given indices, with `inds` a tuple of index ranges or the ranges given separately as `xs` and `ys`, e.g.

```julia
model = gaussian(x=0, y=0, fwhm=(8, 10))
psfplot(model, (-30:30, -30:30))
```

These functions have no methods until a plotting package is loaded:

  - `using Plots` enables a [Plots.jl](https://github.com/JuliaPlots/Plots.jl) user recipe, so all Plots keyword arguments are supported.
  - `using Makie` (e.g. via GLMakie/CairoMakie) enables a [Makie.jl](https://github.com/MakieOrg/Makie.jl) recipe, so all `heatmap` attributes are supported.

Loading both backends in the same session is not supported.
"""
function psfplot end

"""
    psfplot!([plt_or_ax,] model, args...; kwargs...)

Mutating variant of [`psfplot`](@ref), plotting into the current (or given) plot or axis.
"""
function psfplot! end

"""
    psfplotview(model, inds; kwargs...)
    psfplotview(model, xs, ys; kwargs...)
    psfplotview(fig[i, j], model, args...; kwargs...)

A [`psfplot`](@ref) bundled with a colorbar (Makie only, requires Makie 0.25+). Returns a block holding the created axis and plot in the `ax` and `plt` fields, e.g.

```julia
fig, view = psfplotview(model, (-30:30, -30:30); colorscale=log10)
view.ax.title = "PSF"
```

Supported keyword arguments include the [`psfplot`](@ref) attributes plus `colorbar` (show the colorbar, default `true`), `colorbar_label`, and `axis` (attributes forwarded to the created Axis).
"""
function psfplotview end

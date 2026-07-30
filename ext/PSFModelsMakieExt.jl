# Makie plotting extension for PSFModels.jl (requires Makie 0.25+).
#
# Implements `psfplot`/`psfplot!` as a Makie compute-graph recipe: the model
# is evaluated over the given indices and displayed with `heatmap`, so the
# recipe supports all of Makie's colormapping attributes and works natively
# with `Makie.Colorbar`. Also implements `psfplotview`, a complex recipe block
# (new in Makie 0.25) bundling an Axis and a Colorbar, analogous to what the
# Plots.jl recipe displays by default.
module PSFModelsMakieExt

using PSFModels
using Makie
using Makie: Aspect, Auto, DataAspect
import PSFModels: psfplot, psfplot!, psfplotview

# The tuple/vector `inds` forms accepted in place of separate `xs`, `ys`
const IndsLike = Union{Tuple, AbstractVector{<:AbstractVector}}

unpackinds(inds) =
    length(inds) == 2 ? Tuple(inds) :
    throw(ArgumentError("`psfplot` requires two index ranges. Got: $inds"))

# ---------------------------------------------------------------------------
# psfplot recipe
# ---------------------------------------------------------------------------

"""
    psfplot(model, xs, ys; kwargs...)
    psfplot(model, inds; kwargs...)

Plot the PSF `model` as a heatmap evaluated over the given indices, e.g.

```julia
model = gaussian(x=0, y=0, fwhm=(8, 10))
psfplot(model, -30:30, -30:30; colorscale=log10, colorrange=(1e-5, 1))
```

See also [`psfplotview`](@ref) for a version bundling a colorbar.

## Arguments

- `model`: the PSF model (any callable evaluating to the PSF value at `(x, y)`)
- `xs`, `ys`: the index ranges to evaluate the model over, given either
  separately or as a single tuple `inds = (xs, ys)`
"""
@recipe PsfPlot (psf, xs, ys) begin
    "Sets whether colors should be interpolated between pixels."
    interpolate = false
    Makie.mixin_colormap_attributes()...
    Makie.mixin_generic_plot_attributes()...
end

# Unpack the tuple/vector `inds` form: `psfplot(model, (xs, ys))`
Makie.convert_arguments(::Type{<:PsfPlot}, psf, inds::IndsLike) = (psf, unpackinds(inds)...)

function Makie.plot!(p::PsfPlot)
    attr = p.attributes
    map!(attr, [:psf, :xs, :ys], :psfdata) do psf, xs, ys
        return [float(psf(x, y)) for x in xs, y in ys]
    end
    heatmap!(p, attr, p.xs, p.ys, p.psfdata)
    return p
end

# Axis defaults matching the Plots.jl recipe (applied when the Axis is
# created by this plot).
function Makie.preferred_axis_attributes(::Type{Makie.Axis}, ::PsfPlot)
    return (; aspect = DataAspect(), xlabel = "x", ylabel = "y")
end

# ---------------------------------------------------------------------------
# PsfPlotView: a complex recipe block (new in Makie 0.25) bundling an Axis
# with a Colorbar — what the Plots.jl recipe approximated with its default
# colorbar. Usage: `fig, view = psfplotview(model, inds)` or
# `psfplotview(fig[1, 1], model, inds)`.
# ---------------------------------------------------------------------------

# Gap between the heatmap and the colorbar, in px (Makie's layout default of
# 18 is noticeably loose for an image panel).
const COLORBAR_GAP = 8

@Block PsfPlotView (psf, xs, ys) begin
    ax::Makie.Axis
    plt::Makie.Plot
    @attributes begin
        "Sets whether colors should be interpolated between pixels."
        interpolate = false
        "Display a colorbar."
        colorbar = true
        "Colorbar label."
        colorbar_label = ""
        "Attributes forwarded to the created Axis, overriding the defaults, e.g. `axis = (; title = \"PSF\")`."
        axis = (;)
        Makie.mixin_colormap_attributes()...
    end
end

Makie.convert_arguments(::Type{PsfPlotView}, psf, inds::IndsLike) = (psf, unpackinds(inds)...)

function psfplotview(args...; kwargs...)
    res = PsfPlotView(args...; kwargs...)
    # Called without a figure position, `PsfPlotView` creates its own Figure,
    # whose default size has nothing to do with the model extent: an
    # aspect-locked panel then sits in it surrounded by whitespace. Shrink the
    # figure onto its content instead. An explicit `figure = (; size = ...)`
    # wins, and a view placed into someone else's figure never resizes it.
    if res isa Makie.FigureBlock && !(:size in keys(get(kwargs, :figure, (;))))
        Makie.update_state_before_display!(res.figure)
        if all(!isnothing, res.block.layoutobservables.autosize[])
            # An explicitly sized view (fixed axis width/height) reports its
            # own footprint, so the standard shrink-wrap applies directly.
            Makie.resize_to_layout!(res.figure)
        else
            fittofigure!(res.figure, res.block)
        end
    end
    return res
end

# `resize_to_layout!` cannot do this for us: an `Aspect` row or column leaves the
# enclosing layout's size undetermined, so a block containing one reports no
# preferred size to the figure and the figure has nothing to shrink onto. But the
# gap between the figure and the block's cell (padding and protrusions) is fixed,
# so the size the content wants can be recovered from the block's own layout.
function fittofigure!(fig::Makie.Figure, bl::PsfPlotView)
    # `tight_bbox` would read the inner layout's `suggestedbbox`, which a Block
    # never sets — it drives `align_to_bbox!` from its `computedbbox` instead.
    have = bl.layoutobservables.computedbbox[]
    _, cells = Makie.GridLayoutBase.compute_rowcols(bl.layout, have)
    want = (cells.rights[end] - cells.lefts[1], cells.tops[1] - cells.bottoms[end])
    size = Makie.widths(Makie.viewport(fig.scene)[]) .+ want .- Makie.widths(have)
    all(>(0), size) && Makie.resize!(fig, round.(Int, size)...)
    return fig
end

function Makie.initialize_block!(bl::PsfPlotView)
    axattrs = Dict{Symbol, Any}(:aspect => DataAspect(), :xlabel => "x", :ylabel => "y")
    # User-provided axis attributes override the defaults
    for (k, v) in pairs(bl.axis[])
        axattrs[k] = v
    end
    # e.g. `psfplotview(model, inds; axis = (; height = 400))`: derive the
    # width to match
    derivepanelsize!(axattrs, extentratio(bl.xs[], bl.ys[]))
    ax = Makie.Axis(bl[1, 1]; axattrs...)
    plt = psfplot!(
        ax, bl.psf, bl.xs, bl.ys;
        interpolate = bl.interpolate, colormap = bl.colormap,
        colorscale = bl.colorscale, colorrange = bl.colorrange,
        lowclip = bl.lowclip, highclip = bl.highclip,
        nan_color = bl.nan_color, alpha = bl.alpha,
    )
    cb = nothing
    if bl.colorbar[]
        cb = Makie.Colorbar(bl[1, 2], plt; label = bl.colorbar_label)
        Makie.colgap!(bl.layout, COLORBAR_GAP)
    end
    # An explicitly sized axis already has its shape: the Aspect-based cell
    # shaping below would fight it, and Aspect rows/columns make the layout
    # nondeterminable, which is exactly what fixed sizes are used to avoid.
    sized = haskey(axattrs, :width) || haskey(axattrs, :height)
    ratio = cellaspectratio(axattrs, bl.xs[], bl.ys[])
    (sized || isnothing(ratio)) || lockcellaspect!(bl, ax, cb, ratio)
    bl.ax = ax
    bl.plt = plt
    return
end

# Extent of the heatmap in data coordinates: the cells extend half a step
# beyond the first/last centers.
halfstep(v) = length(v) > 1 ? abs(last(v) - first(v)) / (2 * (length(v) - 1)) : 0.5

function extentratio(xs, ys)
    w = abs(last(xs) - first(xs)) + 2 * halfstep(xs)
    h = abs(last(ys) - first(ys)) + 2 * halfstep(ys)
    return w / h
end

# Fixed axis sizes keep the layout fully determined (so `resize_to_layout!`
# works around any arrangement of panels); an `Aspect`-shaped row or column
# does not. When only one of width/height is fixed, derive the other from the
# data extent, so the panel matches the aspect the data would be displayed at.
function derivepanelsize!(axattrs, ratio)
    w = get(axattrs, :width, nothing)
    h = get(axattrs, :height, nothing)
    isnothing(w) && h isa Real && (axattrs[:width] = h * ratio)
    isnothing(h) && w isa Real && (axattrs[:height] = w / ratio)
    return axattrs
end

# The width:height the *axis box* should have, given the attributes it was
# created with, or `nothing` if it is free to fill its cell.
function cellaspectratio(axattrs, xs, ys)
    aspect = get(axattrs, :aspect, nothing)
    aspect isa Real && return Float64(aspect)
    aspect isa DataAspect || return nothing
    # Explicit user limits win over the data extent
    lims = get(axattrs, :limits, nothing)
    lims isa NTuple{4, Any} && (lims = ((lims[1], lims[2]), (lims[3], lims[4])))
    if lims isa Tuple{Any, Any}
        x, y = lims
        if x isa Tuple{Any, Any} && y isa Tuple{Any, Any} && all(!isnothing, (x..., y...))
            return (x[2] - x[1]) / (y[2] - y[1])
        end
    end
    return extentratio(xs, ys)
end

# An `Axis` `aspect` shrinks the axis *inside* its layout cell, so the cell keeps
# whatever shape the layout gives it and the leftover space shows up as a gap
# between the heatmap and the colorbar. Give the cell itself the right shape
# instead — then the axis fills it and the colorbar sits flush against the image.
# (This is the approach in Makie's "Aspect ratios and automatic figure sizes"
# tutorial.) `Aspect` derives one cell dimension from the other, and the derived
# one is unbounded: deriving the width from the height overflows the layout for
# wide panels, deriving the height from the width overflows it for tall ones. So
# pick the direction that fits in the space the block was actually given, and
# revisit it whenever that space changes.
function lockcellaspect!(bl::PsfPlotView, ax::Makie.Axis, cb, ratio::Real)
    layout = bl.layout
    # The block's bbox is already net of protrusions (tick labels, axis labels,
    # colorbar labels), so the only thing between the axis cell and the block's
    # right edge is the colorbar column and the gap before it.
    function reservedwidth()
        isnothing(cb) && return 0.0
        bar = something(cb.layoutobservables.autosize[][1], 0.0f0)
        gap = max(ax.layoutobservables.protrusions[].right, cb.layoutobservables.protrusions[].left)
        return bar + gap + COLORBAR_GAP
    end
    function relayout(bbox)
        w, h = Makie.widths(bbox)
        availw = max(w - reservedwidth(), 1.0)
        availh = max(h, 1.0)
        colsize, rowsize = if availw / availh >= ratio
            Aspect(1, ratio), Auto()   # height binds: derive the width from it
        else
            Auto(), Aspect(1, 1 / ratio)   # width binds: derive the height from it
        end
        (layout.colsizes[1] == colsize && layout.rowsizes[1] == rowsize) && return
        # Both must change together: an Aspect row and an Aspect column cannot be
        # resolved at the same time, so the layout errors on the intermediate state.
        Makie.GridLayoutBase.with_updates_suspended(layout) do
            layout.colsizes[1] = colsize
            layout.rowsizes[1] = rowsize
        end
        return
    end
    relayout(bl.layoutobservables.computedbbox[])
    Makie.on(relayout, bl.layoutobservables.computedbbox)
    return
end

end # module

# Plots.jl extension for PSFModels.jl, implementing `psfplot`/`psfplot!` as a
# RecipesBase user recipe.
module PSFModelsPlotsExt

using PSFModels
using RecipesBase
import PSFModels: psfplot, psfplot!

# The tuple/vector `inds` forms accepted in place of separate `xs`, `ys`
const IndsLike = Union{Tuple, AbstractVector{<:AbstractVector}}

mutable struct PsfPlot
    args
end

# What `RecipesBase.@userplot PsfPlot` would generate, except with concrete
# signatures (instead of `args...`) attached to the stubs in PSFModels, so the
# methods cannot collide with the ones the Makie extension defines.
psfplot(model, inds::IndsLike; kw...) = RecipesBase.plot(PsfPlot((model, inds)); kw...)
psfplot(model, xs::AbstractVector, ys::AbstractVector; kw...) = RecipesBase.plot(PsfPlot((model, xs, ys)); kw...)
psfplot!(model, inds::IndsLike; kw...) = RecipesBase.plot!(PsfPlot((model, inds)); kw...)
psfplot!(model, xs::AbstractVector, ys::AbstractVector; kw...) = RecipesBase.plot!(PsfPlot((model, xs, ys)); kw...)
psfplot!(plt::RecipesBase.AbstractPlot, args...; kw...) = RecipesBase.plot!(plt, PsfPlot(args); kw...)

@recipe function f(p::PsfPlot)
    model = p.args[1]
    inds = p.args[2:end]
    # if `inds` is a vector/tuple, this effectively unpacks it
    if length(inds) == 1
        inds = inds[1]
    end

    seriestype := :heatmap
    aspect_ratio --> 1
    xlims --> extrema(first(inds))
    ylims --> extrema(last(inds))
    xguide --> "x"
    yguide --> "y"

    arr = map(model, CartesianIndices(inds))

    return inds..., transpose(arr)
end

end # module

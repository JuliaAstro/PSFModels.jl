using RecipesBase

@recipe function f(model::PSFModel, inds...)
    seriestype := :heatmap
    aspect_ratio --> 1
    xlims --> extrema(first(inds))
    ylims --> extrema(last(inds))

    arr = model[reverse(inds)...]

    return inds..., arr
end

@recipe f(model::PSFModel, inds=reverse(axes(model))) = model, inds...

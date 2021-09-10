using RecipesBase

@recipe function f(model::PSFModel, inds...)
    seriestype := :heatmap
    aspect_ratio --> 1
    xlims --> extrema(last(inds))
    ylims --> extrema(first(inds))

    arr = model[inds...]

    return reverse(inds)..., arr
end

@recipe f(model::PSFModel, inds=axes(model)) = model, inds...

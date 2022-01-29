using RecipesBase

@recipe function f(model::PSFModel, inds...)
    seriestype := :heatmap
    aspect_ratio --> 1
    xlims --> extrema(first(inds))
    ylims --> extrema(last(inds))
    xguide --> "x"
    yguide --> "y"

    arr = model[inds...]

    return inds..., transpose(arr)
end

@recipe f(model::PSFModel, inds=axes(model)) = model, inds...

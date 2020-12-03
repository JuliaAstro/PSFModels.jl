using RecipesBase

@recipe function f(model::PSFModel, inds...)
    seriestype := :heatmap
    aspect_ratio --> 1
    xlim --> extrema(inds[2])
    ylim --> extrema(inds[1])

    arr = model[inds...]

    return inds..., arr
end

@recipe f(model::PSFModel, inds=axes(model)) = model, inds...

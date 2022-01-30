using RecipesBase

@userplot PsfPlot

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

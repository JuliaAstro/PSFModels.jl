using RecipesBase

@userplot PsfPlot

PsfPlot(model, inds...) = PsfPlot(model, Tuple(inds))

@recipe function f(p::PsfPlot)
    model, inds = p.args

    seriestype := :heatmap
    aspect_ratio --> 1
    xlims --> extrema(first(inds))
    ylims --> extrema(last(inds))
    xguide --> "x"
    yguide --> "y"

    arr = model.(CartesianIndices(inds))

    return inds..., transpose(arr)
end

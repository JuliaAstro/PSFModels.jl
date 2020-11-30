using RecipesBase

@recipe function f(kernel::PSFKernel)
    seriestype := :heatmap
    aspect_ratio --> 1

    inds = axes(kernel)
    arr = kernel[inds...]

    xlim --> extrema(inds[2])
    ylim --> extrema(inds[1])

    return inds..., arr
end

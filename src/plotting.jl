using RecipesBase

@recipe f(kernel::PSFKernel, sz=size(kernel))
    seriestype := :heatmap
    aspect_ratio --> 1

    inds = CartesianIndices(sz)
    arr = kernel[inds]

    xlim --> extrema(axes(arr, 2))
    ylim --> extrema(axes(arr, 1))

    return arr
end

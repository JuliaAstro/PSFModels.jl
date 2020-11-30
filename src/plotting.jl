using RecipesBase

@recipe function f(kernel::PSFKernel, inds...)
    seriestype := :heatmap
    aspect_ratio --> 1
    xlim --> extrema(inds[2])
    ylim --> extrema(inds[1])

    arr = kernel[inds...]

    return inds..., arr
end

@recipe f(kernel::PSFKernel, inds=axes(kernel)) = kernel, inds...

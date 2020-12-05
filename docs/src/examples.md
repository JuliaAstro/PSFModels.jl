
# Examples

## Fitting a PSF

Here is a brief example which shows how to construct a loss function for fitting a `PSFModel` to some data.

```@example fit
using PSFModels: Gaussian
using HCIDatasets: BetaPictoris
using Plots

# convenience function for plotting
function imshow(data; kwargs...)
    xlim = extrema(axes(data, 2))
    ylim = extrema(axes(data, 1))
    heatmap(data; xlim=xlim, ylim=ylim, aspect_ratio=1, kwargs...)
end

# get a PSF from HCIDatasets.jl;
# you may be prompted to download the file
psf = BetaPictoris[:psf]

imshow(psf)
```

```@example fit
using LossFunctions

# generative model
function model(X::AbstractVector{T}) where T
    position = @view X[1:2] # x, y position
    fwhm     = @view X[3:4] # fwhm_x, fwhm_y
    amp      =       X[5]   # amplitude
    return amp * Gaussian{T}(position, fwhm)
end

# objective function
function loss(X::AbstractVector{T}, target) where T
    # cheap way to enforce positivity
    all(>(0), X) || return T(Inf)
    # get generative model
    m = model(X)
    # l2-distance loss (χ² loss) (LossFunctions.jl)
    stamp = @view m[axes(target)...]
    return value(L2DistLoss(), target, stamp, AggMode.Sum())
end

# params are [x, y, fwhm_x, fwhm_y, amp]
test_params = Float32[20, 20, 5, 5, 1]
loss(test_params, psf)
```

The objective function can then be used with an optimization library like [Optim.jl](https://github.com/JuliaOpt/Optim.jl) to find best-fitting parameters

```@example fit
using Optim

# Fit our data using test_params as a starting point
# uses Nelder-Mead optimization
res = optimize(P -> loss(P, psf), test_params)
```

```@example fit
# utilize automatic differentiation (AD) to enable
# advanced algorithms, like LBFGS
res_ad = optimize(P -> loss(P, psf), test_params, LBFGS(); autodiff=:forward)
```

we can see which result has the better loss, and then use the generative model to create a model that we can use elsewhere

```@example fit
best_res = minimum(res) < minimum(res_ad) ? res : res_ad
best_fit_params = Optim.minimizer(best_res)
```

```@example fit
synth_psf = model(best_fit_params)

plot(
    imshow(psf, title="Data"),
    plot(synth_psf, axes(psf); title="Model"),
    cbar=false,
    ticks=false,
    layout=2,
    size=(600, 300)
)
```

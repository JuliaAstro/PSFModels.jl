
# Examples

## Fitting a PSF

Here is a brief example which shows how to construct a loss function for fitting a `PSFModel` to some data.

```@example fit
using PSFModels
using HCIDatasets: BetaPictoris
using Plots
using Statistics

# convenience function for plotting
function imshow(data; kwargs...)
    xlim = extrema(axes(data, 1))
    ylim = extrema(axes(data, 2))
    heatmap(transpose(data); xlim=xlim, ylim=ylim, aspect_ratio=1, kwargs...)
end

# get a PSF from HCIDatasets.jl;
# you may be prompted to download the file
psf = BetaPictoris[:psf]

imshow(psf)
```

```@example fit
# generative model
function model(X::AbstractVector{T}) where T
    x    =       X[1]   # position
    y    =       X[2]
    fwhm = @view X[3:4] # fwhm_x, fwhm_y
    amp  =       X[5]   # amplitude
    return airydisk(T; x, y, fwhm, amp)
end

# objective function
function loss(X::AbstractVector{T}, target) where T
    # cheap way to enforce positivity
    all(>(0), X) || return T(Inf)
    # get generative model
    m = model(X)
    # mean square error
    return mean(idx -> (m(idx) - psf[idx])^2, CartesianIndices(psf))
end

# params are [x, y, fwhm_x, fwhm_y, amp]
test_params = eltype(psf)[20, 20, 5, 5, 1]
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
    psfplot(synth_psf, axes(psf); title="Model"),
    cbar=false,
    ticks=false,
    xlabel="",
    ylabel="",
    layout=2,
    size=(600, 300)
)
```

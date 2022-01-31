
# Examples

## Fitting a PSF

Here is a brief example which shows how to construct a loss function for fitting a `PSFModel` to some data.

```@example fit
using PSFModels
using PSFModels: fit
using HCIDatasets: BetaPictoris
using Plots
using Statistics

# convenience function for plotting
function imshow(data; kwargs...)
    xlim = extrema(axes(data, 1))
    ylim = extrema(axes(data, 2))
    heatmap(transpose(data); xlim=xlim, ylim=ylim,
            aspect_ratio=1, clims=(1e-5, Inf), kwargs...)
end

# get a PSF from HCIDatasets.jl;
# you may be prompted to download the file
psf = BetaPictoris[:psf]

imshow(psf)
```

We can fit this data with a variety of models, here showcasing the flexible [`PSFModels.fit`](@ref) function.


### Gaussian

Using [`gaussian`](@ref)

```@example fit
# parameter vector must match order of values we want to fit
params = (:x, :y, :fwhm, :amp)
P0 = Float32[20, 20, 5, 0.1]
P_gauss, mod_gauss = fit(gaussian, params, P0, psf)
pairs(P_gauss)
```

```@example fit
plot(
    imshow(psf, title="Data"),
    imshow(mod_gauss, title="Model"),
    cbar=false,
    ticks=false,
    xlabel="",
    ylabel="",
    layout=2,
    size=(600, 300)
)
```

and now using a rotated, elliptical Gaussian

```@example fit
# parameter vector must match order of values we want to fit
params = (:x, :y, :fwhm, :amp, :theta)
# use two values for theta here, one for each axis
P0 = Float32[20, 20, 5, 5, 0.1, 0]
P_ellip, mod_ellip = fit(gaussian, params, P0, psf)
pairs(P_ellip)
```

```@example fit
plot(
    imshow(psf, title="Data"),
    imshow(mod_ellip, title="Model"),
    cbar=false,
    ticks=false,
    xlabel="",
    ylabel="",
    layout=2,
    size=(600, 300)
)
```

### Airy disk

Now with [`airydisk`](@ref)

```@example fit
# parameter vector must match order of values we want to fit
params = (:x, :y, :fwhm, :amp, :ratio)
P0 = Float32[20, 20, 5, 0.1, 0.3]
P_airy, mod_airy = fit(airydisk, params, P0, psf)
pairs(P_airy)
```

```@example fit
plot(
    imshow(psf, title="Data"),
    imshow(mod_airy, title="Model"),
    cbar=false,
    ticks=false,
    xlabel="",
    ylabel="",
    layout=2,
    size=(600, 300)
)
```

### Moffat

And finally, with [`moffat`](@ref)


```@example fit
# parameter vector must match order of values we want to fit
params = (:x, :y, :fwhm, :amp, :theta, :alpha)
# again, two values for fwhm for each axis
P0 = Float32[20, 20, 5, 5, 0.1, 0, 2]
P_moff, mod_moff = fit(moffat, params, P0, psf)
pairs(P_moff)
```

```@example fit
plot(
    imshow(psf, title="Data"),
    imshow(mod_moff, title="Model"),
    cbar=false,
    ticks=false,
    xlabel="",
    ylabel="",
    layout=2,
    size=(600, 300)
)
```

### Changing optimization parameters

Any keyword arguments get passed on to `Optim.optimize`, and you can change the algorithm used with the `alg` keyword

```@example fit
using Optim

# parameter vector must match order of values we want to fit
params = (:x, :y, :fwhm, :amp, :theta, :alpha)
# again, two values for fwhm for each axis
P0 = Float32[20, 20, 5, 5, 0.1, 0, 2]
P_moff, mod_moff = fit(moffat, params, P0, psf; alg=Newton())
pairs(P_moff)
```


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
params = (x=20, y=20, fwhm=5, amp=0.1)
P_gauss, mod_gauss = fit(gaussian, params, psf)
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
params = (x=20, y=20, fwhm=(5, 5), amp=0.1, theta=0)
P_ellip, mod_ellip = fit(gaussian, params, psf)
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
params = (x=20, y=20, fwhm=5, amp=0.1, ratio=0.3)
P_airy, mod_airy = fit(airydisk, params, psf)
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
params = (x=20, y=20, fwhm=(5, 5), amp=0.1, theta=0, alpha=2)
P_moff, mod_moff = fit(moffat, params, psf)
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
# load Optim.jl to use the Newton method
using Optim

params = (x=20, y=20, fwhm=(5, 5), amp=0.1, theta=0, alpha=2)
P_moff, mod_moff = fit(moffat, params, psf; alg=Newton())
pairs(P_moff)
```

We can also "freeze" parameters by creating a named tuple and passing it to `func_kwargs`

```@example fit

params = (;x=10, y=20, fwhm=(5, 5), amp=0.1)
func_kwargs = (;alpha=2)
P_moff2, mod_moff2 = fit(moffat, params, psf; func_kwargs)
pairs(P_moff2)
```

```@example fit
plot(
    imshow(psf, title="Data"),
    imshow(mod_moff2, title="Model"),
    cbar=false,
    ticks=false,
    xlabel="",
    ylabel="",
    layout=2,
    size=(600, 300)
)
```
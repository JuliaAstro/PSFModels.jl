# API/Reference

```@setup plots
using PSFModels
using Plots
default(colorbar_scale=:log10, clims=(1e-5, 1))
```

```@index
```

## Gaussian

```@docs
gaussian
normal
```

```@example plots
gauss = gaussian(x=0, y=0, fwhm=10)
psfplot(gauss, -50:50, -50:50; title="gaussian(fwhm=10)")
```

## Airy Disk

```@docs
airydisk
```

```@example plots
airy = airydisk(x=0, y=0, fwhm=10)
psfplot(airy, -50:50, -50:50; title="airydisk(fwhm=10)")
```

## Moffat

```@docs
moffat
```

```@example plots
moff = moffat(x=0, y=0, fwhm=10)
psfplot(moff, -50:50, -50:50; title="moffat(fwhm=10)")
```

## Comparison

```@example plots
xs = range(0, 50, length=100)
plot(
    xs, [gauss.(xs, 0) airy.(xs, 0) moff.(xs, 0)], 
    label=["gaussian" "airydisk" "moffat"], yscale=:log10,
    xlabel="x", ylabel="I/I0", ylims=(1e-5, 1)
)
```
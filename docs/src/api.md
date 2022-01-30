# API/Reference

```@setup plots
using PSFModels
using Plots
```

```@index
```

## Gaussian

```@docs
gaussian
normal
```

```@example plots
model = gaussian(x=0, y=0, fwhm=10)
psfplot(model, (-30:30, -30:30); title="gaussian(fwhm=10)")
```

## Airy Disk

```@docs
airydisk
```

```@example plots
model = airydisk(fwhm=10)
psfplot(model, (-30:30, -30:30); title="airydisk(fwhm=10)")
```

## Moffat

```@docs
moffat
```

```@example plots
model = moffat(fwhm=10)
psfplot(model, (-30:30, -30:30); title="moffat(fwhm=10)")
```

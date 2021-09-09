# API/Reference

```@setup plots
using PSFModels: Gaussian, Moffat, AiryDisk
using Plots
```

```@index
```

```@docs
PSFModels.PSFModel
```

## Gaussian

```@docs
PSFModels.Gaussian
PSFModels.Normal
```

```@example plots
model = Gaussian(10)
plot(model; title="Gaussian(fwhm=10)")
```

## Airy Disk

```@docs
PSFModels.AiryDisk
```

```@example plots
model = AiryDisk(10)
plot(model; title="AiryDisk(fwhm=10)")
```

## Moffat

```@docs
PSFModels.Moffat
```

```@example plots
model = Moffat(10)
plot(model; title="Moffat(fwhm=10)")
```

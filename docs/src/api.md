# API/Reference

```@setup plots
using PSFModels: Gaussian, Moffat, AiryDisk
using Plots
```

```@index
```

```@docs
PSFModels.PSFKernel
```

## Gaussian

```@docs
PSFModels.Gaussian
PSFModels.Normal
```

```@example plots
kernel = Gaussian(10)
plot(kernel; title="Gaussian(fwhm=10)")
```

## Airy Disk

```@docs
PSFModels.AiryDisk
```

```@example plots
kernel = AiryDisk(10)
plot(kernel; title="AiryDisk(fwhm=10)")
```

## Moffat

```@docs
PSFModels.Moffat
```

```@example plots
kernel = Moffat(10)
plot(kernel; title="Moffat(fwhm=10)")
```

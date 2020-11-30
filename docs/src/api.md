# API/Reference

```@setup plots
using PSFKernels: Gaussian, Moffat, AiryDisk
using Plots
```

```@index
```

```@docs
PSFKernels.PSFKernel
```

## Gaussian

```@docs
PSFKernels.Gaussian
PSFKernels.Normal
```

```@example plots
kernel = Gaussian(10)
plot(kernel; title="Gaussian(fwhm=10)")
```

## Airy Disk

```@docs
PSFKernels.AiryDisk
```

```@example plots
kernel = AiryDisk(10)
plot(kernel; title="AiryDisk(fwhm=10)")
```

## Moffat

```@docs
PSFKernels.Moffat
```

```@example plots
kernel = Moffat(10)
plot(kernel; title="Moffat(fwhm=10)")
```

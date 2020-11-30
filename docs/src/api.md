# API/Reference

```@setup plots
using PSFKernels: Gaussian, Moffat, AiryDisk
using Plots
# convenience function for plotting
function imshow(data; kwargs...)
    xlim = extrema(axes(data, 2))
    ylim = extrema(axes(data, 1))
    heatmap(data; xlim=xlim, ylim=ylim, aspect_ratio=1, kwargs...)
end
```

```@index
```

```@docs
PSFKernels.PSFKernel
```

# Gaussian

```@eval
kernel = Gaussian(10)
imshow(kernel; title="Gaussian(10; maxsize=3))
```

```@docs
PSFKernels.Gaussian
PSFKernels.Normal
```

## Airy Disk

```@eval
kernel = AiryDisk(10)
imshow(kernel; title="AiryDisk(10; maxsize=3))
```

```@docs
PSFKernels.AiryDisk
```

## Moffat

```@eval
kernel = Moffat(10)
imshow(kernel; title="Moffat(10; maxsize=3))
```

```@docs
PSFKernels.Moffat
```

```@meta
CurrentModule = PSFModels
```

# Effective PSF Models

An effective point-spread function (ePSF) describes how a point source appears
after both the optical PSF and the detector pixel response have acted on the
incoming light. Following [Anderson2000](@citet), if the instrumental PSF is
$\psi_I$ and the pixel response is $\mathcal{R}$, the ePSF is

```math
\psi_E = \psi_I \ast \mathcal{R}.
```

With this convention, a background-subtracted detector pixel can be modeled as

```math
P_{ij} - s_\ast =
f_\ast\,\psi_E(i - x_\ast, j - y_\ast),
```

where $(x_\ast, y_\ast)$ is the stellar centroid, $f_\ast$ is the stellar flux,
and $s_\ast$ is the local background. The ePSF therefore predicts detector-pixel
fluxes directly from offsets between pixel centers and the source centroid,
without requiring a separate pixel integral at every model evaluation.

## Why Use an ePSF?

The ePSF formulation is especially useful when the PSF is comparable to, or
narrower than, the detector pixel scale. In that regime, a star's measured pixel
values depend strongly on its subpixel location. By combining many stars with
different pixel phases, an empirical ePSF can recover a higher-resolution model
of the pixel-convolved image of a point source.

In `PSFModels.jl`, [`ImagePSF`](@ref) can fit an ePSF from cutouts of different
stars on the same image, each of which will have different pixel phases, allowing
for the measurement of an ePSF. It stores this model as an image sampled on an
oversampled grid and evaluates it with bicubic interpolation. This makes it a
practical model for fitting stellar centroids and fluxes while preserving the
detector-pixel flux convention used by the analytic PSF and PRF models.

## Single-Image Limitations

Single-image ePSF measurement is still suboptimal because the fitted centroids
and the reconstructed ePSF are not independent. The builder alternates between
using the current ePSF to fit each star's centroid and using those fitted
centroids to project the same detector pixels back onto the oversampled ePSF
grid. Any pixel-phase-dependent centroid error in one step can therefore be
absorbed into the next ePSF estimate, leaving an irreducible systematic
degeneracy between the shape of the ePSF and the centroids used to construct it.

This is the limitation addressed by the multi-image treatment in Section 4.4 of
[Anderson2000](@citet). When the same stars are observed in multiple carefully
chosen dithers, each star lands at different pixel phases in different images.
The shared stellar positions then provide an external constraint on the
centroids, while the shifted pixel phases provide complementary samples of the
same underlying ePSF. Considering the dithered images simultaneously breaks much
of the single-image degeneracy and gives a better constrained ePSF than any one
image can provide on its own.

The current [`ImagePSF`](@ref) builder is a single-image implementation, so it
uses the ensemble of stars in one exposure to estimate the ePSF and anchors only
the median centroid drift. Support for simultaneous multi-image, dither-aware
ePSF construction is planned for future development.

## Available Models

- [`ImagePSF`](@ref): an empirical, image-backed ePSF model built from isolated
  stars in a single image. See [ImagePSF](@ref image_psf) for the builder API
  and implementation details.

## Oversampling and number of PSF stars

Oversampling increases the number of ePSF grid cells that must be constrained
from the same detector pixels. A useful first-order estimate comes from treating
the stellar pixel phases as independent uniform random variables on
$[0, 1) \times [0, 1)$. Let the oversampling be
$(s_x, s_y)$ and define

```math
q = s_x s_y
```

as the number of oversampled phase cells per detector pixel. For one ePSF grid
cell away from cutout boundaries, each star contributes one relevant detector
pixel whose random pixel phase lands in that cell with probability

```math
p = \frac{1}{q} = \frac{1}{s_x s_y}.
```

With $N_\ast$ usable PSF stars, the number of detector-pixel samples assigned to
that ePSF cell is therefore

```math
K \sim \mathrm{Binomial}\!\left(N_\ast, \frac{1}{s_x s_y}\right),
```

so

```math
\mathbb{E}[K] = \frac{N_\ast}{s_x s_y}.
```

Equivalently, to obtain an average of $\mu$ input samples per ePSF cell,

```math
N_\ast \approx \mu\,s_x s_y.
```

This scaling is fairly intuitive; doubling the oversampling along both axes
requires roughly four times as many stars for the same per-cell support. It also
gives the expected hole rate before interpolation. For an interior ePSF cell,

```math
\Pr(K = 0) = \left(1 - \frac{1}{s_x s_y}\right)^{N_\ast},
```

and requiring an expected empty-cell fraction below $f$ gives

```math
N_\ast \ge
\frac{\log f}{\log\left(1 - \frac{1}{s_x s_y}\right)}
\approx s_x s_y \log\frac{1}{f}.
```

The empty-cell criterion is weaker than the requirement for a high-quality
stack. For example, a 4x oversampled grid needs only about 72 stars to make the
expected interior hole fraction less than 1%, but that corresponds to only
$72 / 16 = 4.5$ samples per cell on average. That is usually enough for the
hole filler to be mostly inactive, but it is still a sparse sample for robust
median combination and sigma rejection.

As a rule of thumb, use at least $\mu \approx 5$ samples per ePSF cell for
exploratory work and $\mu \approx 10$ or more for production-quality stacks,
assuming isolated stars with good S/N. For scalar oversampling $s = s_x = s_y$:

| Oversampling | Phase cells $s^2$ | Stars for $\mu=5$ | Stars for $\mu=10$ | Stars for <1% holes |
|---:|---:|---:|---:|---:|
| 1× | 1 | 5 | 10 | 1 |
| 2× | 4 | 20 | 40 | 16 |
| 4× | 16 | 80 | 160 | 72 |
| 8× | 64 | 320 | 640 | 293 |
| 16× | 256 | 1280 | 2560 | 1177 |

These counts should be interpreted as usable stars after masking, clipping, and
quality cuts. Boundary cells, bad pixels, saturation masks, cosmic-ray rejection,
crowding, low S/N, and non-uniform pixel phases all reduce the effective number
of samples, so real data often need more stars than this idealized calculation.
If the build reports many filled holes or the fitted ePSF changes noticeably
when the star list is resampled, reduce the oversampling or add more stars.

## References

```@bibliography
Pages = ["epsf_overview.md"]
Canonical = false
```

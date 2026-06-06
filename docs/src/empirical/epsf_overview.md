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
  stars in a single image. See [ImagePSF](image_psf.md) for the builder API,
  implementation details, and guidance on choosing the oversampling factor and
  number of PSF stars.

## References

```@bibliography
Pages = ["epsf_overview.md"]
Canonical = false
```

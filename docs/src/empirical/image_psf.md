```@meta
CurrentModule = PSFModels
```

# ImagePSF — ePSF model for single images

## Overview

[`ImagePSF`](@ref) is an empirical, image-backed effective
point-spread function (PSF) model. 
Unlike the analytic PSF/PRF models ([`GaussianPSF`](@ref), [`AiryPSF`](@ref), etc.),
`ImagePSF` represents the PSF as a discrete 2D array of samples on an
oversampled grid relative to the detector pixel scale. It is evaluated at
arbitrary sub-pixel locations via separable bicubic interpolation.

The distinction between an *instrumental* PSF and an *effective* PSF (ePSF) is
central to the `ImagePSF` design. Following [Anderson2000](@citet), the effective
PSF $\psi_E(\Delta x, \Delta y)$ is the convolution of the instrumental PSF
$\psi_I$ with the pixel response function $\mathcal{R}$, so that a
background-subtracted stellar pixel value can be written directly as

```math
P_{ij} - s_\ast = f_\ast \, \psi_E(i - x_\ast,\, j - y_\ast)
```

without performing a pixel integral at every evaluation. This is the convention
used throughout `ImagePSF`: the tabulated `data` grid samples $\psi_E$, and the
model's `flux` parameter directly scales detector-pixel fluxes.

---

## Public API

### Building an ePSF from data

```@docs
fit(::Type{ImagePSF}, ::AbstractMatrix, ::Any, ::Any; fit_rad)
fit(::Type{ImagePSF}, ::AbstractMatrix, ::Any; x, y)
```

### Result types

```@docs
ImagePSF
```

```@docs
ImagePSFBuildResult
```

---

## Supporting Methods

The following methods are part of the `AbstractPSFModel` interface and work
identically for `ImagePSF` as for analytic models:

| Method | Description |
|---|---|
| `evaluate(model, x, y)` | Evaluate the ePSF at detector-pixel coordinates |
| `centroid(model)` | Return `(model.x, model.y)` |
| `integral(model)` | Return `model.flux` |
| `background(model)` | Return `model.bkg` |
| `extent(model)` | Bounding box in detector-pixel coordinates |
| `effective_area(model)` | Reciprocal of the sum of squared ePSF samples, scaled by oversampling area |

---

## Implementation

The ePSF construction algorithm in [`build_epsf`](@ref) (called by
`fit(ImagePSF, ...)`) is an adaptation of the [`Anderson2000`](@citet) iterative
residual-stacking method to the **single-image** regime. The classical method
operates across dithered multi-exposure datasets with three stages: (1) ePSF
construction, (2) per-star fitting, and (3) multi-dither position averaging to
break the PSF/centroid degeneracy. This implementation covers (1) and (2),
with (3) planned in the future.

### Algorithm outline

The builder iterates between two phases — stacking the ePSF grid and fitting
individual stars — until the largest centroid shift falls below `centroid_tol`:

```
1. EXTRACT  → extract_stars(image, x, y, fit_rad)
2. STACK    → stack_epsf_grid(image, stars, state)
3. FIT      → fit_all_stars(stars, psf, image)
4. ANCHOR   → remove_centroid_drift(stars, old_centroids)
5. RESTACK  → stack_epsf_grid(image, stars, state)
6. Repeat 3–5 until convergence
```

### Step 1 — Star extraction

[`extract_stars`](@ref) takes an image and initial detector-pixel coordinates
`(x, y)`, creates a square cutout of radius `fit_rad` around each position, and
keeps only pixels whose full square area lies inside that radius. For each star
an [`EmpiricalStar`](@ref) record is created with:

- The cutout index ranges and a flat vector of pixel coordinates
- An initial sky estimate from [`estimate_local_sky`](@ref) (median of the
  cutout's edge pixels)
- An initial flux estimate from [`estimate_initial_star_params`](@ref), which
  uses a winsorized (98th-percentile-capped) aperture sum of
  background-subtracted pixels to suppress hot defects

**Comparison to [Anderson2000](@citet):** Anderson & King use a modal sky from
an annulus (4–7 pixel radii) and a more sophisticated source-detection pipeline
(DAOPHOT or aperture photometry). Our simpler window-edge sky estimate works for
the isolated-star regime but would be biased in crowded fields. The winsorized
flux is a pragmatic defense against cosmic rays and hot pixels that the
classical method handles through its sigma-rejection stacking step.
**Plan to integrate more closely with Photometry.jl to improve this step.**

```@docs
PSFModels.extract_stars
PSFModels.EmpiricalStar
PSFModels.estimate_local_sky
PSFModels.estimate_initial_star_params
```

### Step 2 — ePSF grid stacking

[`stack_epsf_grid`](@ref) builds one iteration's ePSF grid from the current
star parameters. It composes six sub-steps:

**2a. Pixel projection:** [`project_star_pixels_to_grid`](@ref) converts each
valid, finite detector pixel into a normalized ePSF sample

```math
\text{sample} = \frac{P_{ij} - s_\ast}{f_\ast}
```

and maps it to the oversampled grid via

```math
g_x = \mathrm{round}\bigl(\mathrm{origin}_x + \mathrm{os}_x \cdot (i - x_\ast)\bigr)
```

using "round half away from zero" to avoid half-grid bias. Samples with
$|\text{sample}| >$ `sample_clip` are rejected as likely pixel defects.

**2b. Robust combination:** [`robust_combine_grid_cells`](@ref) applies a
sigma-clipped median — the median is computed, then samples beyond
`sigma_clip` × 1.4826 × MAD are rejected and the surviving median is kept.
Cells with no samples are left as `NaN`.

```@docs
PSFModels.robust_combine_grid_cells
```

**Comparison to [Anderson2000](@citet):** The classical method forms
*residuals* (sample minus current ePSF model), averages them within a 0.25-pixel
radius of each grid point, rejects at 2.5σ, and adjusts the grid point by
the mean residual, iterating this inner loop 5 times. Our method directly
stacks normalized samples rather than residuals, uses a median statistic
(rather than mean), and does not iterate the stacking inner loop. The
direct-stack approach requires fewer stacking iterations but depends on
having a good initial flux estimate. The median is more robust to outliers
than the sigma-clipped mean, which may be advantageous when fewer training
stars are available.

**2c. Hole filling:** [`fill_grid_holes!`](@ref) repairs `NaN` cells in
three stages: (i) interior holes with finite support on both sides along both
axes are filled by separable cubic Lagrange interpolation, iterating up to
`maxiter=6` passes by default so that newly filled cells can support their neighbors;
(ii) border holes that lack bracketing support are extrapolated with
Anderson's local exponential radial-tail fit (a linear fit of
$\log(\mathrm{value})$ vs. radius from the grid center, using positive
samples in a 5×5 neighborhood); (iii) any remaining holes are filled with
local neighbor medians, falling back to zero.

```@docs
PSFModels.fill_grid_holes!
```

**Comparison to [Anderson2000](@citet):** The bicubic infill stage has no 
direct counterpart in [`Anderson2000`](@citet) and may introduce mild
smoothing artifacts in very sparsely sampled regions, but if more than 10%
of ePSF pixels are empty, a warning is issued. In well-sampled grids it is
effectively a no-op.

**2d. Smoothing:** [`smooth_grid_quartic!`](@ref) convolves the grid with
the fixed 5×5 quartic kernel of [Anderson2000](@citet) (their Eq. 8).
This step is identical to the classical method when `smooth=true`.

**2e. Recentering:** [`recenter_grid_to_origin!`](@ref) estimates the
grid's mass-weighted centroid and shifts the grid back to the nominal
origin via bicubic interpolation. Shifts larger than 2 oversampling pixels
are ignored as likely stack failures.

```@docs
PSFModels.smooth_grid_quartic!
PSFModels.recenter_grid_to_origin!
```

**Comparison to [Anderson2000](@citet):** The classical method uses an
edge-symmetry condition: it adjusts the origin so that
$\psi_E(-0.5, \delta y) = \psi_E(+0.5, \delta y)$ and similarly for $y$,
using derivatives at the central-pixel edges. Our mass-weighted centroid is
simpler but can be biased by asymmetric grid coverage (e.g., when most stars
sample one side of the grid more densely). For well-distributed stellar
samples the two should agree closely.

**2f. Normalization:** [`normalize_grid_to_oversampling_area!`](@ref)
optionally clips negative values and rescales the grid so that
$\sum \mathrm{data} = \mathrm{os}_x \cdot \mathrm{os}_y$. This enforces the
convention that a unit-flux star's ePSF integrates to one detector-pixel
flux per oversampling cell — the same convention used by
[`evaluate`](@ref) for [`ImagePSF`](@ref).

```@docs
PSFModels.normalize_grid_to_oversampling_area!
```

### Step 3 — Per-star fitting

[`fit_star_against_epsf`](@ref) fits each star's centroid $(x, y)$, flux, and
optionally background against the current ePSF. It uses Levenberg-Marquardt
optimization ([`lm_irls`](@ref)) with analytic derivatives provided by
[`evaluate_fg`](@ref). When `reweight` is a robust loss (e.g.,
[`TukeyLoss`](@ref)), the fit uses iteratively reweighted least squares to
down-weight outlier pixels.

```@docs
PSFModels.fit_star_against_epsf
```

**Comparison to [Anderson2000](@citet):** The classical method uses a
Newton-Raphson solver with analytic derivatives derived from a quadratic
expansion of $\chi^2$, operating on a soft-edged aperture of ~1.5 pixels
with Poisson variance weighting. Our LM/IRLS approach is more general
(supports arbitrary loss functions and optional background fitting), uses
the full star cutout rather than a small aperture, and does not incorporate
Poisson variance weighting by default (though `inv_var` can be passed to
`fit_lm` if needed). The larger fitting aperture makes our method less
sensitive to the exact choice of `fit_rad` but somewhat more computationally
expensive per star.

### Step 4 — Centroid anchoring

[`remove_centroid_drift`](@ref) computes the median centroid shift across all
surviving stars and subtracts it. This breaks the global degeneracy between
the PSF shape and the absolute centroid reference frame.

```@docs
PSFModels.remove_centroid_drift
```

**Comparison to [Anderson2000](@citet):** This is the single most significant
deviation from the classical method. [Anderson2000](@citet) Stage 3 uses
multiple dithered exposures: the same star observed at different pixel phases
provides an external constraint that pins down the absolute astrometric frame.
Our median-drift removal is a single-image heuristic — it prevents the PSF
from wandering but cannot correct pixel-phase error that is common to all
stars. In dithered datasets, the classical multi-exposure averaging is
unambiguously superior. In truly single-image contexts, median anchoring is
the best one can do without an external reference. The method will be most
reliable when the training stars uniformly sample pixel phase; systematic
centroid biases can remain if stars cluster at particular phases.
**We plan to add support for analyzing multipel dithered exposures
together in the future**.

### Convergence

Convergence is achieved when the largest accepted
centroid update across all stars falls below `centroid_tol` (default
$10^{-3}$ pixels). A minimum of two surviving stars is required.

### Return value

[`build_epsf`](@ref) returns a `(psf, result)` pair, where `result` is an
[`ImagePSFBuildResult`](@ref) containing the final per-star parameters,
usage and convergence flags, iteration count, and per-star costs — enabling
downstream quality filtering.

## References
This page cites the following references:

```@bibliography
Pages = ["image_psf.md"]
Canonical = false
```
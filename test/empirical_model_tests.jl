using PSFModels
using PSFModels: _as_oversampling, background, bicubic_interpolate, centroid, evaluate_fg, extent, fit_lm, integral, fill_grid_holes!, TukeyLoss
import ConstructionBase
using StableRNGs: StableRNG
using Statistics: mean
using Test

"""
    _truth_grid(model, psf::ImagePSF)

Evaluate an analytic unit-flux model on the same oversampled grid as an
`ImagePSF`, then normalize it to the `ImagePSF` sum convention.
"""
function _truth_grid(model, psf::ImagePSF)
    # Convert each oversampled grid point into detector-pixel offsets.
    truth = similar(psf.data)
    ox, oy = psf.origin
    sx, sy = psf.oversampling
    for j in axes(truth, 2), i in axes(truth, 1)
        truth[i, j] = evaluate(model, (i - ox) / sx, (j - oy) / sy)
    end
    # Match ImagePSF's normalization for direct array comparisons.
    truth .*= prod(psf.oversampling) / sum(truth)
    return truth
end

@testset "ImagePSF API and bicubic interpolation" begin
    # Check core model API behavior because ImagePSF stores fixed grid data
    # alongside fit parameters, unlike the analytic models.
    data = [exp(-((i - 4)^2 + (j - 4)^2) / 5) for i in 1:7, j in 1:7]
    model = ImagePSF(data; x = 10.3, y = 11.4, flux = 120.0, bkg = 7.0, oversampling = 2, normalize = true)

    @test centroid(model) == (10.3, 11.4)
    @test integral(model) == 120.0
    @test background(model) == 7.0
    @test sum(model.data) ≈ 4.0
    @test extent(model) == ((8.8, 11.8), (9.9, 12.9))
    @test evaluate(model, -100, -100) == 7.0

    fill_model = ImagePSF(data; x = 0, y = 0, flux = 4, bkg = 1, fill_value = 0.25)
    @test evaluate(fill_model, -100, -100) == 2.0

    # Updating fit parameters should reuse the immutable model's PSF array
    # instead of making a copy.
    updated = ConstructionBase.setproperties(model, (x = 1, flux = 2))
    @test updated.x === 1.0
    @test updated.flux === 2.0
    @test updated.data === model.data # Identity check for array reuse

    f, g = evaluate_fg(model, 10.8, 11.9)
    @test f ≈ evaluate(model, 10.8, 11.9)

    # Compare the analytic ImagePSF gradient against finite differences.
    p0 = [model.x, model.y, model.flux, model.bkg]
    fd = similar(p0)
    h = 1.0e-5
    for k in eachindex(p0)
        pplus = copy(p0)
        pminus = copy(p0)
        pplus[k] += h
        pminus[k] -= h
        mplus = ImagePSF(model.data; x = pplus[1], y = pplus[2], flux = pplus[3], bkg = pplus[4], origin = model.origin, oversampling = model.oversampling)
        mminus = ImagePSF(model.data; x = pminus[1], y = pminus[2], flux = pminus[3], bkg = pminus[4], origin = model.origin, oversampling = model.oversampling)
        fd[k] = (evaluate(mplus, 10.8, 11.9) - evaluate(mminus, 10.8, 11.9)) / (2h)
    end
    @test collect(g) ≈ fd rtol = 1.0e-5 atol = 1.0e-5
end

@testset "bicubic_interpolate" begin
    data = [exp(-((i - 4)^2 + (j - 4)^2) / 5) for i in 1:7, j in 1:7]
    model = ImagePSF(data; x = 10.3, y = 11.4, flux = 120.0, bkg = 7.0, oversampling = 2, normalize = true)
    # Verify out-of-bounds interpolation and constructor validation paths.
    val, dx, dy = bicubic_interpolate(model.data, -1, 2; fill_value = 0.3)
    @test val == 0.3
    @test dx == 0
    @test dy == 0
    @test_throws ArgumentError bicubic_interpolate(rand(3, 4), 2, 2)
    @test_throws ArgumentError ImagePSF(rand(3, 4))
    @test_throws ArgumentError ImagePSF(data; oversampling = (2.0, 2))
    bad = copy(data)
    bad[1, 1] = NaN
    @test_throws ArgumentError ImagePSF(bad)
end

@testset "fill_grid_holes!" begin
    # Infill behavior
    x=rand(StableRNG(43), 21,21)
    # inds=(vcat(rand(1:21, 5), 1), vcat(rand(1:21, 5), 1))
    inds = ([9, 4, 15, 13, 1, 1], [6, 15, 17, 1, 17, 1])
    x[inds...] .= NaN
    fill_grid_holes!(x)
    @test all(isfinite, x)
    # Regression: bordered by three good pixels, should not copy single neighbor
    @test x[1,1] != x[1,2]
    # Test that warning triggers if too much data is missing,
    # but still fills what it can without error.
    fill!(x, 1)
    inds = (1:2:21, 1:2:21)
    x[inds...] .= NaN
    @test_logs (:warn,) fill_grid_holes!(x)
    @test all(==(1), x) # Fallback infill should fill all holes in this case
end

@testset "_as_oversampling" begin
    # Confirm accepted scalar and axis-specific oversampling forms.
    @test _as_oversampling(2) == (2, 2)
    @test _as_oversampling((2, 3)) == (2, 3)
    @test _as_oversampling([3, 4]) == (3, 4)

    # Reject non-integer, non-positive, and incorrectly sized inputs.
    @test_throws ArgumentError _as_oversampling(2.0)
    @test_throws ArgumentError _as_oversampling((2.0, 2))
    @test_throws ArgumentError _as_oversampling((2, 0))
    @test_throws ArgumentError _as_oversampling((2, 3, 4))
end

@testset "ImagePSF LM fit -- single star" begin
    # Fit only ImagePSF's source parameters to ensure it works with LM/IRLS.
    grid_model = CircularGaussianPRF(x = 8, y = 8, fwhm = 2.4, flux = 1, bkg = 0)
    psf_data = render(grid_model, (1:16, 1:16))
    truth = ImagePSF(psf_data; x = 8.35, y = 7.75, flux = 300.0, bkg = 4.0, origin = (8.0, 8.0), normalize = true)
    image = render(truth, (1:16, 1:16))
    init = ImagePSF(psf_data; x = 8.0, y = 8.1, flux = 260.0, bkg = 3.5, origin = (8.0, 8.0), normalize = true)

    best, result = fit_lm(init, image, (1:16, 1:16); max_iter = 100)
    @test result.converged
    @test best.x ≈ truth.x atol = 1e-6
    @test best.y ≈ truth.y atol = 1e-6
    @test best.flux ≈ truth.flux rtol = 1e-6
    @test best.bkg ≈ truth.bkg atol = 1e-6
end

@testset "simulation utilities" begin
    # Exercise source placement, deterministic rendering, SNR conversion, and noise.
    rng = StableRNG(10)
    psf = CircularGaussianPRF(x = 0, y = 0, fwhm = 2.5, flux = 1, bkg = 0)
    sources = simulate_sources((40, 45), 8; flux = (100.0, 120.0), min_separation = 4, border = 5, rng)
    @test length(sources.x) == 8
    @test all(6 .≤ sources.x .≤ 35)
    @test all(6 .≤ sources.y .≤ 40)
    for i in 1:length(sources.x), j in (i + 1):length(sources.x)
        @test hypot(sources.x[i] - sources.x[j], sources.y[i] - sources.y[j]) ≥ 4
    end

    image = simulate_image((40, 45), psf, sources; background = 12.0, noise = :none, model_radius = 7)
    @test size(image) == (40, 45)
    @test sum(image) ≈ 12.0 * 40 * 45 + sum(sources.flux) rtol = 1.0e-6
    @test flux_for_snr(psf; snr = 20, background = 12, read_noise = 2) > 0

    noisy = copy(image)
    add_noise!(noisy; noise = :poisson_gaussian, read_noise = 1.0, rng)
    @test noisy != image
end

@testset "Empirical ImagePSF recovery" begin
    # Recover an undersampled integrated Gaussian PRF from many clean stars.
    rng = StableRNG(42)
    truth_model = CircularGaussianPRF(x = 0, y = 0, fwhm = 1.8, flux = 1, bkg = 0)
    image, sources = simulate_image(
        (128, 128),
        truth_model,
        70;
        background = 20.0,
        noise = :none,
        flux = (600.0, 900.0),
        min_separation = 10,
        border = 8,
        model_radius = 6,
        rng
    )

    psf, result = PSFModels.fit(
        ImagePSF,
        image,
        sources.x .+ randn(rng, length(sources.x)) * 0.5,
        sources.y .+ randn(rng, length(sources.y)) * 0.5;
        fit_rad = 5.0,
        oversampling = 2,
        smooth = true,
        recenter = true,
        anchor_centroids = false,
        reweight = nothing,
    )
    truth = _truth_grid(truth_model, psf)
    @test sum(psf.data) ≈ 4.0
    @test count(result.used) ≥ 65
    @test mean(abs.(psf.data .- truth)) < 0.003
    @test mean(hypot.(sources.x[result.used] .- result.x[result.used], sources.y[result.used] .- result.y[result.used])) < 0.05
    # Guard against extrapolating unconstrained edge holes into a false PSF tail.
    @test maximum(abs.(psf.data .- truth)[20:21, 15:21]) < 1e-5

    # Fit the stars in the image with the measured PSF model to confirm it can reproduce the input fluxes and centroids.
    for k in eachindex(sources.x)
        # model = ImagePSF(psf.data; x = sources.x[k], y = sources.y[k], flux = sources.flux[k], origin = psf.origin, oversampling = psf.oversampling)
        model = ConstructionBase.setproperties(psf, (x = sources.x[k], y = sources.y[k], flux = sources.flux[k]))
        fit, _ = fit_lm(model, image, (1:128, 1:128))
        @test fit.flux ≈ sources.flux[k] rtol = 0.05
        @test fit.x ≈ sources.x[k] atol = 0.1
        @test fit.y ≈ sources.y[k] atol = 0.1
    end

    # Verify the explicit-cutout API forwards into the same empirical builder.
    inds = ntuple(
        k -> (
            floor(Int, sources.x[k] - 5):ceil(Int, sources.x[k] + 5),
            floor(Int, sources.y[k] - 5):ceil(Int, sources.y[k] + 5),
        ),
        length(sources.x)
    )
    psf_from_inds, inds_result = PSFModels.fit(
        ImagePSF,
        image,
        inds;
        x = sources.x,
        y = sources.y,
        oversampling = 2,
        maxiters = 1,
        smooth = true,
        recenter = true,
        anchor_centroids = false,
        reweight = nothing,
    )
    @test count(inds_result.used) == 70
    @test sum(psf_from_inds.data) ≈ 4.0
end

@testset "Empirical ImagePSF with defective stars" begin
    # Stress the robust stack with many hot and low pixels in the training stars.
    rng = StableRNG(123)
    truth_model = CircularGaussianPRF(x = 0, y = 0, fwhm = 1.9, flux = 1, bkg = 0)
    image, sources = simulate_image(
        (96, 96),
        truth_model,
        70;
        background = 15.0,
        noise = :none,
        flux = (600.0, 900.0),
        min_separation = 7,
        border = 8,
        model_radius = 6,
        rng
    )

    # Contaminate a high fraction of stars to test defect rejection.
    for k in 1:round(Int, 0.75 * length(sources.x))
        dx = rand(rng, -3:3)
        dy = rand(rng, -3:3)
        i = clamp(round(Int, sources.x[k]) + dx, 1, size(image, 1))
        j = clamp(round(Int, sources.y[k]) + dy, 1, size(image, 2))
        image[i, j] += rand(rng) < 0.8 ? 1500.0 : -14.0
    end

    psf, result = PSFModels.fit(
        ImagePSF,
        image,
        sources.x .+ randn(rng, length(sources.x)) * 0.25,
        sources.y .+ randn(rng, length(sources.y)) * 0.25;
        fit_rad = 5.0,
        oversampling = 2,
        smooth = true,
        recenter = true,
        anchor_centroids = false,
    )
    truth = _truth_grid(truth_model, psf)
    @test count(result.used) ≥ 65
    @test mean(abs.(psf.data .- truth)) < 1e-3
    @test mean(hypot.(sources.x[result.used] .- result.x[result.used], sources.y[result.used] .- result.y[result.used])) < 0.075
end

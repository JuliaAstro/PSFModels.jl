using PSFModels
using StableRNGs: StableRNG
using Test

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

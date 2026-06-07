using PSFModels
using Distributions: Poisson
using StableRNGs: StableRNG
using Test

"""
    ks_2sample_statistic(x, y)

Return the two-sample Kolmogorov–Smirnov statistic between samples `x` and `y`.

This computes

    D = sup_z |F_x(z) - F_y(z)|

where `F_x` and `F_y` are empirical CDFs.
"""
function ks_2sample_statistic(x::AbstractVector, y::AbstractVector)
    n = length(x)
    m = length(y)
    n > 0 || throw(ArgumentError("x must be nonempty"))
    m > 0 || throw(ArgumentError("y must be nonempty"))

    xs = sort(x)
    ys = sort(y)

    i = j = 1
    cdf_x = 0.0
    cdf_y = 0.0
    d = 0.0

    while i <= n || j <= m
        if j > m || (i <= n && xs[i] < ys[j])
            z = xs[i]
        elseif i > n || ys[j] < xs[i]
            z = ys[j]
        else
            z = xs[i]
        end

        while i <= n && xs[i] == z
            i += 1
        end
        cdf_x = (i - 1) / n

        while j <= m && ys[j] == z
            j += 1
        end
        cdf_y = (j - 1) / m

        d = max(d, abs(cdf_x - cdf_y))
    end

    return d
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

    @testset "flux_for_snr" begin
        # Verify that requested SNR inverts the data-unit noise model, including source shot noise.
        snr = 20.0
        background = 12.0
        read_noise = 2.0
        gain = 2.5
        flux = flux_for_snr(psf; snr, background, read_noise, gain)
        σb2 = background / gain + read_noise^2
        @test flux / sqrt(flux / gain + effective_area(psf) * σb2) ≈ snr

        # Source-only Poisson noise still requires nonzero flux in data units.
        @test flux_for_snr(psf; snr = 5, background = 0, gain = 2) ≈ 12.5
    end

    noisy = copy(image)
    add_noise!(noisy; noise = :poisson_gaussian, read_noise = 1.0, rng)
    @test noisy != image

    @testset "poisson" begin
        # Compare the internal Poisson sampler against Distributions.jl across small and large rates.
        for val in (0.5, 2.0, 10.0, 1000.0)
            rng = StableRNG(123)
            x = [PSFModels.rand_poisson(rng, val) for _ in 1:100_000]
            y = rand(rng, Poisson(val), 100_000)
            @test ks_2sample_statistic(x, y) < 0.01
        end
    end
end

using PSFModels
using PSFModels: free_params, model_from_vector, _has_hessian, _has_deriv, fit, fit_lm
using PSFModels: MADScale, FixedScale, MScale, TukeyLoss, estimate_scale, weight
using Distributions: Poisson
import LossFunctions
import Optim
using LinearAlgebra: diag
using StableRNGs
using Statistics: mean, median, std
using Test

# ---------------------------------------------------------------------------
# Tests for free_params / model_from_vector
# ---------------------------------------------------------------------------

const irls_losses = (
    LossFunctions.HuberLoss(1.0),
    TukeyLoss(),
    nothing, # no IRLS ⇒ no reweighting
)

@testset "free_params / model_from_vector" begin
    m = CircularGaussianPSF(x=1.0, y=2.0, fwhm=4.0, flux=10.0, bkg=1.0)
    names, idx, x0 = free_params(m)
    @test names == (:x, :y, :fwhm, :flux, :bkg)
    @test idx == (1, 2, 3, 4, 5)
    @test x0 == [1.0, 2.0, 4.0, 10.0, 1.0]

    names_fixed, idx_fixed, x0_fixed = free_params(m, (bkg=1.0,))
    @test names_fixed == (:x, :y, :fwhm, :flux)
    @test idx_fixed == (1, 2, 3, 4)
    @test length(x0_fixed) == 4

    m2 = model_from_vector(m, Val(names), [1.1, 2.2, 3.3, 9.0, 0.5], (;))
    @test m2.x ≈ 1.1
    @test m2.y ≈ 2.2
    @test m2.fwhm ≈ 3.3
    @test m2.flux ≈ 9.0
    @test m2.bkg ≈ 0.5
end

# ---------------------------------------------------------------------------
# Tests for capability detection
# ---------------------------------------------------------------------------

@testset "_has_hessian / _has_deriv" begin
    m1 = CircularGaussianPSF(x=0.0, y=0.0, fwhm=3.0, flux=1.0, bkg=0.0)
    m2 = GaussianPSF(x=0.0, y=0.0, x_fwhm=4.0, y_fwhm=3.0, theta=0.0, flux=1.0, bkg=0.0)
    @test _has_hessian(m1)
    @test _has_deriv(m1)
    @test _has_hessian(m2)
    @test _has_deriv(m2)
end

# ---------------------------------------------------------------------------
# Noiseless recovery tests for CircularGaussianPSF
# ---------------------------------------------------------------------------

@testset "struct fit CircularGaussianPSF" begin
    inds = (1:30, 1:30)
    truth = CircularGaussianPSF(x=13.5, y=12.3, fwhm=4.0, flux=200.0, bkg=5.0)
    img = render(truth, inds)
    # Initial guess: perturbed ~5% away from truth
    init = CircularGaussianPSF(x=14.1, y=11.8, fwhm=4.3, flux=190.0, bkg=5.5)

    # Default L2 / Newton path
    best, result = fit(init, img, inds; x_abstol=1e-6)
    @test best.x    ≈ truth.x    rtol=1e-3
    @test best.y    ≈ truth.y    rtol=1e-3
    @test best.fwhm ≈ truth.fwhm rtol=1e-3
    @test best.flux ≈ truth.flux rtol=1e-3
    @test best.bkg  ≈ truth.bkg  rtol=1e-3
    @test Optim.converged(result)

    # With inv_var — also returns covariance
    inv_var = fill(1.0, size(img))
    best2, result2, cov = fit(init, img, inds; inv_var, x_abstol=1e-6)
    @test best2.x ≈ truth.x rtol=1e-3
    @test size(cov) == (5, 5)

    # Freeze background: init with wrong bkg, but bkg is fixed so it must converge anyway
    init_fixed = CircularGaussianPSF(x=14.1, y=11.8, fwhm=4.3, flux=190.0, bkg=5.0)
    best3, _ = fit(init_fixed, img, inds; fixed=(bkg=5.0,), x_abstol=1e-6)
    @test best3.bkg ≈ 5.0
    @test best3.x ≈ truth.x rtol=1e-3

    # LogCoshLoss (still twice differentiable → Newton path)
    best4, _ = fit(init, img, inds; loss=LossFunctions.LogCoshLoss(), x_abstol=1e-6)
    @test best4.x ≈ truth.x rtol=1e-3

    # HuberLoss (not twice differentiable → LBFGS path)
    best5, _ = fit(init, img, inds; loss=LossFunctions.HuberLoss(1.0), x_abstol=1e-6)
    @test best5.x ≈ truth.x rtol=1e-3
end

# ---------------------------------------------------------------------------
# Noiseless recovery tests for GaussianPSF
# ---------------------------------------------------------------------------

@testset "struct fit GaussianPSF" begin
    inds = (1:40, 1:40)
    truth = GaussianPSF(x=18.5, y=17.3, x_fwhm=5.0, y_fwhm=3.5, theta=15.0, flux=300.0, bkg=2.0)
    img = render(truth, inds)
    # Initial guess: perturbed ~5% away from truth
    init = GaussianPSF(x=19.3, y=16.7, x_fwhm=5.3, y_fwhm=3.3, theta=13.0, flux=285.0, bkg=2.2)

    best, result = fit(init, img, inds; x_abstol=1e-6)
    @test best.x      ≈ truth.x      rtol=1e-3
    @test best.y      ≈ truth.y      rtol=1e-3
    @test best.x_fwhm ≈ truth.x_fwhm rtol=1e-3
    @test best.y_fwhm ≈ truth.y_fwhm rtol=1e-3
    @test best.flux   ≈ truth.flux   rtol=1e-3
    @test Optim.converged(result)

    # Freeze theta: init with slightly wrong theta but it's fixed
    init_fixed = GaussianPSF(x=19.3, y=16.7, x_fwhm=5.3, y_fwhm=3.3, theta=15.0, flux=285.0, bkg=2.2)
    best2, _ = fit(init_fixed, img, inds; fixed=(theta=15.0,), x_abstol=1e-6)
    @test best2.theta ≈ 15.0
    @test best2.x ≈ truth.x rtol=1e-3

    # With inv_var
    inv_var = fill(1.0, size(img))
    best3, result3, cov = fit(init, img, inds; inv_var, x_abstol=1e-6)
    @test best3.x ≈ truth.x rtol=1e-3
    @test size(cov) == (7, 7)
end

# ---------------------------------------------------------------------------
# ArgumentError checks
# ---------------------------------------------------------------------------

@testset "struct fit argument errors" begin
    m = CircularGaussianPSF(x=5.5, y=5.2, fwhm=3.0, flux=100.0, bkg=1.0)
    img = render(m, (1:10, 1:10))
    init = CircularGaussianPSF(x=5.8, y=4.9, fwhm=3.2, flux=95.0, bkg=1.1)
    @test_throws ArgumentError fit(init, img; inv_var=zeros(size(img)))
    @test_throws ArgumentError fit(init, img; inv_var=fill(-1.0, size(img)))
    @test_throws ArgumentError fit(init, img; inv_var=fill(1.0, (11, 10)))
end

# ---------------------------------------------------------------------------
# LM fitting (fit_lm) tests
# ---------------------------------------------------------------------------

@testset "fit_lm noiseless recovery" begin
    inds = (1:30, 1:30)
    truth = CircularGaussianPSF(x=13.5, y=12.3, fwhm=4.0, flux=200.0, bkg=5.0)
    img = render(truth, inds)
    init = CircularGaussianPSF(x=14.1, y=11.8, fwhm=4.3, flux=190.0, bkg=5.5)

    # No inv_var
    best, result = fit_lm(init, img, inds)
    @test result.converged
    @test best.x    ≈ truth.x    rtol=1e-3
    @test best.y    ≈ truth.y    rtol=1e-3
    @test best.fwhm ≈ truth.fwhm rtol=1e-3
    @test best.flux ≈ truth.flux rtol=1e-3
    @test size(result.cov) == (5, 5)
    @test isnan(result.σ_final)  # no IRLS ⇒ no scale estimate

    # With unity inv_var
    inv_var = fill(1.0, size(img))
    best2, result2 = fit_lm(init, img, inds; inv_var)
    @test result2.converged
    @test result.minimizer ≈ result2.minimizer
end

@testset "fit_lm IRLS — noiseless parity" begin
    # On clean data, IRLS should recover the same parameters as plain L2
    inds = (1:30, 1:30)
    truth = CircularGaussianPSF(x=13.5, y=12.3, fwhm=4.0, flux=200.0, bkg=5.0)
    img = render(truth, inds)
    init = CircularGaussianPSF(x=14.1, y=11.8, fwhm=4.3, flux=190.0, bkg=5.5)

    best_l2, result_l2 = fit_lm(init, img, inds)
    best_huber, result_huber = fit_lm(init, img, inds; reweight=LossFunctions.HuberLoss(1.0))
    best_tukey, result_tukey = fit_lm(init, img, inds; reweight=TukeyLoss())

    @test result_huber.converged
    @test result_tukey.converged
    @test best_huber.x ≈ truth.x rtol=1e-2
    @test best_huber.y ≈ truth.y rtol=1e-2
    @test best_tukey.x ≈ truth.x rtol=5e-2
    @test best_tukey.y ≈ truth.y rtol=5e-2
end

@testset "fit_lm IRLS — outlier rejection" begin
    inds = (1:30, 1:30)
    truth = CircularGaussianPSF(x=15.0, y=15.0, fwhm=4.0, flux=200.0, bkg=1.0)
    img = render(truth, inds)
    # Add Gaussian noise and scattered severe outliers in the wings
    rng = StableRNG(42)
    img .= img .+ 0.5 .* randn(rng, size(img))
    # Several 20σ outliers: severe enough to bias flux and FWHM
    img[5, 5]   += 10.0
    img[25, 6]  += 10.0
    img[8, 24]  += 10.0
    img[15, 28] += 10.0
    img[22, 12] += 10.0

    init = CircularGaussianPSF(x=15.5, y=14.5, fwhm=3.5, flux=180.0, bkg=1.5)

    best_l2, result_l2 = fit_lm(init, img, inds; max_iter=500)
    @test result_l2.converged

    best_huber, result_huber = fit_lm(init, img, inds; reweight=LossFunctions.HuberLoss(1.0), max_iter=500)
    @test result_huber.converged

    best_tukey, result_tukey = fit_lm(init, img, inds; reweight=TukeyLoss(), max_iter=500)
    @test result_tukey.converged

    # IRLS improves flux accuracy: outliers inflate the wings, which L2
    # compensates for by reducing the total flux; robust weights suppress
    # those outliers and recover a more accurate flux estimate.
    err_l2_flux = abs(best_l2.flux - truth.flux)
    @test abs(best_huber.flux - truth.flux) < err_l2_flux
    @test abs(best_tukey.flux - truth.flux) < err_l2_flux

    # IRLS improves FWHM accuracy: outlier pixels in the wings cause L2 to
    # over-widen the fitted profile; downweighting them tightens the estimate.
    err_l2_fwhm = abs(best_l2.fwhm - truth.fwhm)
    @test abs(best_huber.fwhm - truth.fwhm) < err_l2_fwhm
    @test abs(best_tukey.fwhm - truth.fwhm) < err_l2_fwhm
end

@testset "fit_lm IRLS — Poisson Errors" begin
    @testset "No inverse variance information" begin
        inds = (1:20, 1:20)
        # for flux in (1000.0, 10000.0)
        bkg = 100.0
        for snr in (50.0, 100.0)
            flux = snr * sqrt(bkg)
            v = [10.0, 10.0, 4.0, flux, bkg]
            truth = CircularGaussianPSF(x=v[1], y=v[2], fwhm=v[3], flux=v[4], bkg=v[5])
            img = render(truth, inds)
            img_noisy = similar(img) # buffer to hold Poisson noise realizations
            # Add Poisson noise
            rng = StableRNG(42)
            init = CircularGaussianPSF(x=v[1]+0.5, y=v[2]+0.5, fwhm=v[3]-0.5, flux=0.8*flux, bkg=v[5]*0.9)

            for loss in (LossFunctions.HuberLoss(1.0), TukeyLoss(), nothing)
                N_tests = 1000
                minimizers = Matrix{Float64}(undef, 5, N_tests)
                for i in 1:N_tests
                    img_noisy .= rand.(Ref(rng), Poisson.(img))
                    best, result = fit_lm(init, img_noisy; reweight=loss)
                    @test result.converged
                    minimizers[:, i] = result.minimizer
                end
                minimizer_means = mean(minimizers, dims=2)[:, 1]
                minimizer_medians = median(minimizers, dims=2)[:, 1]
                minimizer_stds = std(minimizers, dims=2)[:, 1]
                # median over many noise realizations should be close to truth
                @test minimizer_medians ≈ v rtol=1e-2
                # fraction of runs where minimizer is within 1σ of truth should be ~68% for each parameter
                within_1σ = sum(abs.(minimizers .- v) .< minimizer_stds, dims=2)[:, 1] ./ N_tests
                # within_1σ is slightly high (~0.69 -- 0.76); the cov calculation when reweighting
                # is somewhat conservative and may overestimate uncertainties
                @test all(within_1σ .>= 0.68 - 0.1) && all(within_1σ .<= 0.68 + 0.1)
            end
        end
    end
end

@testset "TukeyLoss properties" begin
    loss = TukeyLoss(; c=4.685)
    @test loss(0.0, 0.0) == 0.0
    # For small residuals, should behave like L2
    @test LossFunctions.deriv(loss, 0.1, 0.0) ≈ 0.1 atol=1e-3
    # For residuals beyond c, derivative should be zero
    @test LossFunctions.deriv(loss, 10.0, 0.0) ≈ 0.0 atol=1e-6
    # Weight function
    @test weight(loss, 0.0) ≈ 1.0 atol=1e-6
    @test weight(loss, 10.0) == 0.0
    @test weight(loss, 2.0) > 0.0  # within threshold, positive weight
end

@testset "IRLS weight function" begin
    @test weight(LossFunctions.L2DistLoss(), 1.0) == 1.0
    @test weight(LossFunctions.L2DistLoss(), 0.0) ≈ 1.0 atol=1e-6
    @test weight(LossFunctions.HuberLoss(1.0), 0.5) == 1.0
    @test weight(LossFunctions.HuberLoss(1.0), 2.0) < 1.0
end

@testset "fit_lm argument errors" begin
    m = CircularGaussianPSF(x=5.5, y=5.2, fwhm=3.0, flux=100.0, bkg=1.0)
    img = render(m, (1:10, 1:10))
    init = CircularGaussianPSF(x=5.8, y=4.9, fwhm=3.2, flux=95.0, bkg=1.1)
    @test_throws ArgumentError fit_lm(init, img; inv_var=zeros(size(img)))
    @test_throws ArgumentError fit_lm(init, img; inv_var=fill(-1.0, size(img)))
    @test_throws ArgumentError fit_lm(init, img; inv_var=fill(1.0, (11, 10)))
    @test_throws ArgumentError fit_lm(init, img; fixed=(x=1.0, y=2.0, fwhm=3.0, flux=100.0, bkg=1.0))
end

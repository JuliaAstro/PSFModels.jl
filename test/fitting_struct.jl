using PSFModels
using PSFModels: free_params, model_from_vector, _has_hessian, _has_deriv, fit
import LossFunctions
import Optim
using LinearAlgebra: diag
using Test

# ---------------------------------------------------------------------------
# Tests for free_params / model_from_vector
# ---------------------------------------------------------------------------

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

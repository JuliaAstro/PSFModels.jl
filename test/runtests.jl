using Distributions
using PSFModels
using StableRNGs
using StaticArrays
using Test

rng = StableRNG(112358)

modname(::typeof(gaussian)) = "gaussian"
modname(::typeof(airydisk)) = "airydisk"
modname(::typeof(moffat)) = "moffat"

function test_model_interface(K)
    # test defaults
    m = K(x=0, y=0, fwhm=10)
    val_dir = @inferred K(0, 0; x=0, y=0, fwhm=10)
    val_call = @inferred m(0, 0)
    val_tup = @inferred m((0, 0))
    val_vec = @inferred m(SA[0, 0])
    @test val_dir ≈ val_call ≈ val_tup ≈ val_vec ≈ 1

    # test new position
    m = K(x=12, y=13, fwhm=10)
    val_dir = @inferred K(12, 13; x=12, y=13, fwhm=10)
    @test m(12, 13) ≈ val_dir ≈ 1

    # test position off pixel grid
    m = K(x=12.5, y=13.5, fwhm=10)
    val_dir = @inferred K(12.5, 13.5; x=12.5, y=13.5, fwhm=10)
    @test m(12.5, 13.5) ≈ val_dir ≈ 1


    # test diagonal fwhm
    m = K(x=0, y=0, fwhm=(10, 9))
    @test m(3, 5) ≈ K(3, 5; x=0, y=0, fwhm=(10, 9)) ≈ K(3, 5; x=0, y=0, fwhm=SA[10, 9])

    # test rotation
    m = K(x=0, y=0, fwhm=(3, 5))
    m90 = K(x=0, y=0, fwhm=(5, 3), theta=90)
    @test m(7, 8) ≈ m90(7, 8)

    mwarn = K(x=0, y=0, fwhm=10, theta=30)
    name = modname(K)
    expected_log = (:warn, "isotropic $name is not affected by non-zero rotation angle 30")
    @test_logs expected_log mwarn(0, 0)

    # test different amplitude
    m = K(x=0, y=0, amp=2, fwhm=9)
    val_dir = @inferred K(0, 0; x=0, y=0, amp=2.0, fwhm=9)
    @test m(0, 0) ≈ val_dir ≈ 2.0

    # test background
    m = K(x=0, y=0, fwhm=9, bkg=10)
    val_dir = @inferred K(0, 0; x=0, y=0, fwhm=9, bkg=10)
    @test m(0, 0) ≈ val_dir ≈ 11

    # test different type
    m = K(BigFloat; x=0, y=0, fwhm=10)
    val_dir = @inferred K(BigFloat, 0, 0; x=0, y=0, fwhm=10)
    @test val_dir isa BigFloat
    @test m(0, 0) ≈ val_dir ≈ BigFloat(1)
end

@testset "Model Interface - $K" for K in (gaussian, airydisk, moffat)
    test_model_interface(K)
end

@testset "gaussian" begin
    m = gaussian(x=0, y=0, fwhm=10)
    expected = exp(-4 * log(2) * sum(abs2, SA[1, 2]) / 100)
    @test m(1, 2) ≈ expected

    m = gaussian(x=0, y=0, fwhm=(10, 9))
    wdist = (1/10)^2 + (2/9)^2
    expected = exp(-4 * log(2) * wdist)
    @test m(1, 2) ≈ expected

    # test Normal alias
    @test normal(0, 0; x=0, y=0, fwhm=10) === gaussian(0, 0; x=0, y=0, fwhm=10)
end

@testset "CircularGaussianPSF" begin
    m = CircularGaussianPSF(x=0, y=0, fwhm=10, flux=1, bkg=0)
    @test centroid(m) == (0.0, 0.0)
    @test integral(m) == 1.0
    r1 = evaluate(m, 1, 2)
    @test r1 isa Float64
    @test r1 ≈ 0.0076829778398427705
    m = CircularGaussianPSF(x=0, y=0, fwhm=10, flux=1, bkg=10)
    @test evaluate(m, 1, 2) ≈ 0.0076829778398427705 + 10
    @test fit_deriv(m, 1, 2) ≈ [0.0004260347542393244, 0.0008520695084786488, -0.001323578190848892, 0.0076829778398427705, 1.0]
    @test fit_hessian(m, 1, 2) ≈ [-0.0004024103711416015 4.724876619544591e-5 -0.00015860171014686828 0.0004260347542393244 0.0; 4.724876619544591e-5 -0.00033153722184843266 -0.00031720342029373656 0.0008520695084786488 0.0; -0.00015860171014686828 -0.00031720342029373656 0.00031777260218123345 -0.0013235781908488918 0.0; 0.0004260347542393244 0.0008520695084786488 -0.0013235781908488918 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
end

@testset "GaussianPSF" begin
    m = GaussianPSF(x=0, y=0, x_fwhm=10, y_fwhm=6, theta=30, flux=1, bkg=0)
    @test centroid(m) == (0.0, 0.0)
    @test integral(m) == 1.0
    r1 = evaluate(m, 1, 2)
    @test r1 isa Float64
    @test r1 ≈ 0.011881854589938992
    m = GaussianPSF(x=0, y=0, x_fwhm=10, y_fwhm=6, theta=30, flux=1, bkg=10)
    @test evaluate(m, 1, 2) ≈ 0.011881854589938992 + 10
    @test fit_deriv(m, 1, 2) ≈ [-6.269560630626032e-5, 0.0025675279951917654, -0.0009587636050457796, -0.0015172854575571534, 4.7000306666382166e-5, 0.01188185458993899, 1.0]
    @test fit_hessian(m, 1, 2) ≈ [-0.0009513701779307769 0.0004936505237653131 -0.00020789110903162168 0.0003838214631391412 -2.986906138973509e-6 -6.26956063062603e-5 0.0; 0.0004936505237653131 -0.0009825507698117963 -0.0003301242561598406 -0.0009787987399725283 -3.547464474559758e-5 0.0025675279951917654 0.0; -0.00020789110903162168 -0.0003301242561598406 0.0001273559772475686 0.00012243190355172497 1.4950134657622455e-6 -0.0009587636050457795 0.0; 0.0003838214631391412 -0.0009787987399725283 0.00012243190355172497 0.0002922935533913823 -3.048115727009258e-5 -0.001517285457557153 0.0; -2.986906138973509e-6 -3.547464474559758e-5 1.4950134657622455e-6 -3.048115727009258e-5 -5.148866586397458e-7 4.700030666638217e-5 0.0; -6.26956063062603e-5 0.0025675279951917654 -0.0009587636050457795 -0.001517285457557153 4.700030666638217e-5 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    # equal fwhm + theta=0 reduces to CircularGaussianPSF
    mc = CircularGaussianPSF(x=1.5, y=2.5, fwhm=8, flux=3, bkg=0)
    mg = GaussianPSF(x=1.5, y=2.5, x_fwhm=8, y_fwhm=8, theta=0, flux=3, bkg=0)
    @test evaluate(mc, 3, 4) ≈ evaluate(mg, 3, 4)
end

@testset "airydisk" begin
    fwhm = 10
    m = airydisk(x=0, y=0, fwhm=fwhm)
    ld = fwhm / 1.028
    radius = 1.22 * ld
    # first radius is 0
    @test m(radius, 0) ≈ 0 atol=1e-5
    @test m(-radius, 0) ≈ 0 atol=1e-5
    @test m(0, radius) ≈ 0 atol=1e-5
    @test m(0, -radius) ≈ 0 atol=1e-5

    # second radius is 0
    radius2 = 2.23 * ld
    @test m(radius2, 0) ≈ 0 atol=1e-5
    @test m(-radius2, 0) ≈ 0 atol=1e-5
    @test m(0, radius2) ≈ 0 atol=1e-5
    @test m(0, -radius2) ≈ 0 atol=1e-5

    fwhm = (10, 6)
    m = airydisk(x=0, y=0, fwhm=fwhm)
    r1 = fwhm[1] * 1.18677
    r2 = fwhm[2] * 1.18677
    # first radius is 0
    @test m(r1, 0) ≈ 0 atol=1e-5
    @test m(-r1, 0) ≈ 0 atol=1e-5
    @test m(0, r2) ≈ 0 atol=1e-5
    @test m(0, -r2) ≈ 0 atol=1e-5

    # test with ratio
    mratio = airydisk(x=0, y=0, fwhm=10, amp=1, ratio=sqrt(0.5))
    # test attenuation
    @test mratio(0, 0) ≈ 1
    @test mratio(radius, 0) > m(radius, 0)
    @test mratio(-radius, 0) > m(-radius, 0)
    @test mratio(0, radius) > m(0, radius)
    @test mratio(0, -radius) > m(0, -radius)


    mratio = airydisk(x=0, y=0, fwhm=(10, 6), ratio=sqrt(0.5), amp=4)
    r1 = fwhm[1] * 1.18677
    r2 = fwhm[2] * 1.18677
    # test attenuation
    @test mratio(0, 0) ≈ 4
    @test mratio(r1, 0) > m(r1, 0)
    @test mratio(-r1, 0) > m(-r1, 0)
    @test mratio(0, r2) > m(0, r2)
    @test mratio(0, -r2) > m(0, -r2)

    # https://github.com/JuliaAstro/PSFModels.jl/issues/14
    mratio = airydisk(x = 40.5, y=40.5, fwhm=10, ratio=0.2,amp=1)
    @test mratio(40.5, 40.5) ≈ 1
end

@testset "moffat" begin
    m = moffat(x=0, y=0, fwhm=10)
     expected = inv(1 + sum(abs2, SA[1, 2]) / 25)
    @test m(1, 2) ≈ expected

    m = moffat(x=0, y=0, fwhm=(10, 9))
    wdist = (1/5)^2 + (2/4.5)^2
    expected = inv(1 + wdist)
    @test m(1, 2) ≈ expected

    # different alpha
    m = moffat(x=0, y=0, fwhm=10, alpha=2)
    expected = inv(1 + sum(abs2, SA[1, 2] ./  PSFModels._moffat_fwhm_to_gamma(10, 2)))^2
    @test m(1, 2) ≈ expected

    fwhms = randn(rng, 100) .+ 10
    for alpha in [0.5, 1.0, 1.5]
        gammas = PSFModels._moffat_fwhm_to_gamma.(fwhms, alpha)
        fwhm_calc = PSFModels._moffat_gamma_to_fwhm.(gammas, alpha)
        @test fwhms ≈ fwhm_calc
    end
end

include("plotting.jl")
include("fitting.jl")

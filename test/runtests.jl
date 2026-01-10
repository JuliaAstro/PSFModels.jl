using ChainRulesCore
using ChainRulesTestUtils
using Distributions
using FiniteDifferences
using PSFModels
using StableRNGs
using StaticArrays
using Test

ChainRulesCore.debug_mode() = true

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

    test_model_interface(gaussian)

    @testset "gradients" begin
        FiniteDifferences.to_vec(x::Integer) = Bool[], _ -> x
        # have to make sure PSFs are all floating point so tangents don't have type issues
        psf_iso = gaussian(fwhm=10.0, pos=zeros(2))
        psf_tang = Tangent{typeof(gaussian)}(fwhm=rand(rng), pos=rand(rng, 2), amp=rand(rng), indices=NoTangent())
        point = Float64[1, 2]
        test_frule(psf_iso ⊢ psf_tang, point)
        test_rrule(psf_iso ⊢ psf_tang, point)

        psf_diag = gaussian(fwhm=Float64[10, 8], pos=zeros(2))
        psf_tang = Tangent{typeof(gaussian)}(fwhm=rand(rng, 2), pos=rand(rng, 2), amp=rand(rng), indices=NoTangent())
        test_frule(psf_diag ⊢ psf_tang, point)
        test_rrule(psf_diag ⊢ psf_tang, point)
    end
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

    test_model_interface(airydisk)
end

@testset "moffat" begin
    test_model_interface(moffat)

    m = moffat(x=0, y=0, fwhm=10)
    expected = inv(1 + sum(abs2, SA[1, 2]) / 25)
    @test m(1, 2) ≈ expected

    m = moffat(x=0, y=0, fwhm=(10, 9))
    wdist = (1/5)^2 + (2/4.5)^2
    expected = inv(1 + wdist)
    @test m(1, 2) ≈ expected

    # different alpha
    m = moffat(x=0, y=0, fwhm=10, alpha=2)
    expected = inv(1 + sum(abs2, SA[1, 2] ./ PSFModels._moffat_fwhm_to_gamma(10, 2)))^2
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

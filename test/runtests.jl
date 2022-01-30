using PSFModels
using PSFModels: gaussian, normal, airydisk, moffat
using StaticArrays
using Test

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

    # test diagonal fwhm
    m = K(x=0, y=0, fwhm=(10, 9))
    @test m(3, 5) ≈ K(3, 5; x=0, y=0, fwhm=(10, 9)) ≈ K(3, 5; x=0, y=0, fwhm=SA[10, 9])

    # test rotation
    m = K(x=0, y=0, fwhm=(3, 5))
    m90 = K(x=0, y=0, fwhm=(5, 3), theta=90)
    @test m(7, 8) ≈ m90(7, 8)

    # test different amplitude
    m = K(x=0, y=0, amp=2, fwhm=9)
    val_dir = @inferred K(0, 0; x=0, y=0, amp=2.0, fwhm=9)
    @test m(0, 0) ≈ val_dir ≈ 2.0

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


@testset "airydisk" begin
    fwhm = 10
    m = airydisk(x=0, y=0, fwhm=fwhm)
    radius = fwhm * 1.18677
    # first radius is 0
    @test m(radius, 0) ≈ 0 atol=eps(Float64)
    @test m(-radius, 0) ≈ 0 atol=eps(Float64)
    @test m(0, radius) ≈ 0 atol=eps(Float64)
    @test m(0, -radius) ≈ 0 atol=eps(Float64)

    fwhm = (10, 9)
    m = airydisk(x=0, y=0, fwhm=fwhm)
    r1 = fwhm[1] * 1.18677
    r2 = fwhm[2] * 1.18677
    # first radius is 0
    @test m(r1, 0) ≈ 0 atol=eps(Float64)
    @test m(-r1, 0) ≈ 0 atol=eps(Float64)
    @test m(0, r2) ≈ 0 atol=eps(Float64)
    @test m(0, -r2) ≈ 0 atol=eps(Float64)
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
    expected = inv(1 + sum(abs2, SA[1, 2]) / 25)^2
    @test m(1, 2) ≈ expected
end

# include("plotting.jl")

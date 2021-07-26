using PSFModels
using PSFModels: Gaussian, Normal, AiryDisk, Moffat, ScaledPSFModel
using StaticArrays
using Test

function test_model_interface(K)
    # test defaults
    m = K(10)
    @test size(m) == (31, 31)
    @test axes(m) == (-15:15, -15:15)
    @test m.pos ≈ SA[0, 0]
    @test eltype(m) == Float64

    @test m[0, 0] ≈ m(0, 0) ≈ m(SA[0, 0]) ≈ 1
    @test_throws ErrorException m(CartesianIndex(0, 0))
    @test m[-100, -10] ≈ m(-10, -100)
    @test m(1.0, 1.0) ≈ m[1.0, 1.0]
    @test maximum(m) ≈ 1
    @test 0 ≤ minimum(m) ≤ 1
    
    @test m .* ones(10, 10) ≈ m[1:10, 1:10]

    # test new position
    m = K(12, 13, 10)
    @test m == K((12, 13), 10) == K(SA[12, 13], 10)
    @test size(m) == (31, 31)
    @test axes(m) == (-2:28, -3:27)
    @test m.pos ≈ SA[12, 13]
    @test eltype(m) == Float64

    # test diagonal fwhm
    m = K((10, 9))
    @test m == K((10, 9)) == K(SA[10, 9])
    @test size(m) == (29, 31)
    @test axes(m) == (-14:14, -15:15)
    @test m.pos ≈ SA[0, 0]
    @test eltype(m) == Float64

    # test different maxsize
    m = K((10, 9); maxsize=2)
    @test size(m) == (19, 21)
    @test axes(m) == (-9:9, -10:10)
    @test m.pos ≈ SA[0, 0]
    @test eltype(m) == Float64

    # test different type
    m = K{BigFloat}(10)
    @test eltype(m) == BigFloat
    @test m[0, 0] isa BigFloat

    # test scaling
    m = 20 * K(10)
    @test m isa ScaledPSFModel
    @test m == K(10) * 20 == ScaledPSFModel(20, K(10))
    @test m(0, 0) ≈ 20
    m2 = K(10) / 20
    @test m2 isa ScaledPSFModel
    @test m2 == ScaledPSFModel(1 / 20, K(10))
    @test m2(0, 0) ≈ 1 / 20
    # composite
    comp = 20 * m2
    @test comp.amp ≈ 1 # amplitudes are combined; not just function chaining
    @test comp.model === m2.model
    comp2 = comp / 100
    @test comp2.amp ≈ 0.01
    @test comp2.model === m2.model # model has propagated through many
end

@testset "Model Interface - $K" for K in (Gaussian, AiryDisk, Moffat)
    test_model_interface(K)
end

@testset "Gaussian" begin
    m = Gaussian(10)
    expected = exp(-4 * log(2) * sum(abs2, SA[1, 2]) / 100)
    @test m[2, 1] ≈ m(1, 2) ≈ expected

    m = Gaussian((10, 9))
    wdist = (1/10)^2 + (2/9)^2
    expected = exp(-4 * log(2) * wdist)
    @test m[2, 1] ≈ m(1, 2) ≈ expected

    # test Normal alias
    @test Normal(10) === Gaussian(10)
end


@testset "AiryDisk" begin
    m = AiryDisk(10)
    radius = m.fwhm * 1.18677
    # first radius is 0
    @test m(radius, 0) ≈ 0 atol=eps(Float64)
    @test m(-radius, 0) ≈ 0 atol=eps(Float64)
    @test m(0, radius) ≈ 0 atol=eps(Float64)
    @test m(0, -radius) ≈ 0 atol=eps(Float64)

    m = AiryDisk((10, 9))
    r1 = m.fwhm[1] * 1.18677
    r2 = m.fwhm[2] * 1.18677
    # first radius is 0
    @test m(r1, 0) ≈ 0 atol=eps(Float64)
    @test m(-r1, 0) ≈ 0 atol=eps(Float64)
    @test m(0, r2) ≈ 0 atol=eps(Float64)
    @test m(0, -r2) ≈ 0 atol=eps(Float64)
end

@testset "Moffat" begin
    m = Moffat(10)
    expected = inv(1 + sum(abs2, SA[1, 2]) / 25)
    @test m[2, 1] ≈ m(1, 2) ≈ expected

    m = Moffat((10, 9))
    wdist = (1/5)^2 + (2/4.5)^2
    expected = inv(1 + wdist)
    @test m[2, 1] ≈ m(1, 2) ≈ expected
end

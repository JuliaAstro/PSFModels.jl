using PSFModels
using PSFModels: Gaussian, Normal, AiryDisk, Moffat
using StaticArrays
using Test

function test_model_interface(K)
    # test defaults
    k = K(10)
    @test size(k) == (31, 31)
    @test axes(k) == (-15:15, -15:15)
    @test k.pos ≈ SA[0, 0]
    @test eltype(k) == Float64

    @test k[0, 0] ≈ k(0, 0) ≈ k(SA[0, 0]) ≈ k(CartesianIndex(0, 0)) ≈ 1
    @test k[-100, -100] ≈ k(-100, -100)
    @test_throws ArgumentError k[1.0, 1.0]
    @test maximum(k) ≈ 1
    @test 0 ≤ minimum(k) ≤ 1
    
    @test k .* ones(10, 10) ≈ k[1:10, 1:10]

    # test new position
    k = K(12, 13, 10)
    @test k == K((12, 13), 10) == K(SA[12, 13], 10)
    @test size(k) == (31, 31)
    @test axes(k) == (-3:27, -2:28)
    @test k.pos ≈ SA[12, 13]
    @test eltype(k) == Float64

    # test diagonal fwhm
    k = K((10, 9))
    @test k == K((10, 9)) == K(SA[10, 9])
    @test size(k) == (31, 29)
    @test axes(k) == (-15:15, -14:14)
    @test k.pos ≈ SA[0, 0]
    @test eltype(k) == Float64

    # test different maxsize
    k = K((10, 9); maxsize=2)
    @test size(k) == (21, 19)
    @test axes(k) == (-10:10, -9:9)
    @test k.pos ≈ SA[0, 0]
    @test eltype(k) == Float64

    # test different type
    k = K{BigFloat}(10)
    @test eltype(k) == BigFloat
    @test k[0, 0] isa BigFloat

end

@testset "Model Interface - $K" for K in (Gaussian, AiryDisk, Moffat)
    test_model_interface(K)
end

@testset "Gaussian" begin
    k = Gaussian(10)
    expected = exp(-4 * log(2) * sum(abs2, SA[1, 2]) / 100)
    @test k[1, 2] ≈ expected

    k = Gaussian((10, 9))
    wdist = (1/10)^2 + (2/9)^2
    expected = exp(-4 * log(2) * wdist)
    @test k[1, 2] ≈ expected

    # test Normal alias
    @test Normal(10) === Gaussian(10)
end


@testset "AiryDisk" begin
    k = AiryDisk(10)
    radius = k.fwhm * 1.18677
    # first radius is 0
    @test k(radius, 0) ≈ 0 atol=eps(Float64)
    @test k(-radius, 0) ≈ 0 atol=eps(Float64)
    @test k(0, radius) ≈ 0 atol=eps(Float64)
    @test k(0, -radius) ≈ 0 atol=eps(Float64)

    k = AiryDisk((10, 9))
    r1 = k.fwhm[1] * 1.18677
    r2 = k.fwhm[2] * 1.18677
    # first radius is 0
    @test k(r1, 0) ≈ 0 atol=eps(Float64)
    @test k(-r1, 0) ≈ 0 atol=eps(Float64)
    @test k(0, r2) ≈ 0 atol=eps(Float64)
    @test k(0, -r2) ≈ 0 atol=eps(Float64)
end

@testset "Moffat" begin
    k = Moffat(10)
    expected = inv(1 + sum(abs2, SA[1, 2]) / 25)
    @test k[1, 2] ≈ expected

    k = Moffat((10, 9))
    wdist = (1/5)^2 + (2/4.5)^2
    expected = inv(1 + wdist)
    @test k[1, 2] ≈ expected
end
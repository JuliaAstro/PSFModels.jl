using PSFKernels
using PSFKernels: Gaussian, Normal, AiryDisk, Moffat
using StaticArrays
using Test

function test_kernel_interface(K)
    # test defaults
    k = K(10)
    @test size(k) == (31, 31)
    @test axes(k) == (-15:15, -15:15)
    @test k.pos ≈ SA[0, 0]
    @test eltype(k) == Float64

    @test k[0, 0] ≈ k(0, 0) ≈ k(SA[0, 0]) ≈ 1
    @test k[-100, -100] ≈ k(-100, -100)
    @test_throws ArgumentError k[1.0, 1.0]
    @test maximum(k) ≈ 1
    @test -1 ≤ minimum(k) ≤ 1

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

@testset "Kernel Interface - $K" for K in (Gaussian, AiryDisk, Moffat)
    test_kernel_interface(K)
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
    
end

@testset "Moffat" begin
    
end
using ChainRulesCore
using ChainRulesTestUtils
using FiniteDifferences
using PSFModels
using PSFModels: Gaussian, Normal, AiryDisk, Moffat
using StableRNGs
using StaticArrays
using Test

ChainRulesCore.debug_mode() = true

rng = StableRNG(122)

function test_model_interface(K)
    # test defaults
    m = @inferred K(fwhm=10)
    @test size(m) == (31, 31)
    @test axes(m) == (-15:15, -15:15)
    @test m.pos ≈ SA[0, 0]
    @test eltype(m) == Float64

    val_idx = @inferred m[0, 0]
    val_call = @inferred m(0, 0)
    val_tup = @inferred m(SA[0, 0])
    @test val_idx ≈ val_call ≈ val_tup ≈ 1
    @test_throws ErrorException m(CartesianIndex(0, 0))
    # out of bounds but it's cool
    @test m[-100, -10] ≈ m(-10, -100)
    @test m(1.0, 1.0) ≈ m[1, 1]
    @test maximum(m) ≈ 1
    @test 0 ≤ minimum(m) ≤ 1
    
    # broadcasting to an array that's smaller
    @test m .* ones(10, 10) ≈ ones(10, 10) .* m ≈ m[1:10, 1:10]
    # broadcasting to a larger array
    @test m .* ones(30, 30) ≈ ones(30, 30) .* m ≈ m[1:30, 1:30]

    # test new position
    m = K(x=12, y=13, fwhm=10)
    @test m == K(pos=(12, 13), fwhm=10) == K(pos=SA[12, 13], fwhm=10) == K(pos=[12, 13], fwhm=10)
    @test size(m) == (31, 31)
    @test axes(m) == (-2:28, -3:27)
    @test m.pos ≈ SA[12, 13]
    @test eltype(m) == Float64

    # test polar position
    m = K(r=2, theta=90, fwhm=10)
    @test m == K(ρ=2, θ=90, fwhm=10)
    @test size(m) == (31, 31)
    @test m.pos ≈ SA[0, 2]

    # test diagonal fwhm
    m = K(fwhm=(10, 9))
    @test m == K(fwhm=(10, 9)) == K(fwhm=SA[10, 9])
    @test size(m) == (29, 31)
    @test axes(m) == (-14:14, -15:15)
    @test m.pos ≈ SA[0, 0]
    @test eltype(m) == Float64

    # test different maxsize
    m = K(fwhm=(10, 9); maxsize=2)
    @test size(m) == (19, 21)
    @test axes(m) == (-9:9, -10:10)
    @test m.pos ≈ SA[0, 0]
    @test eltype(m) == Float64

    # test different amplitude
    m = K(amp=2, fwhm=9)
    @test m == K(amp=2.0, fwhm=9)
    @test m.amp == 2.0
    @test m(m.pos) == 2.0

    # test different type
    m = K(BigFloat, fwhm=10)
    @test eltype(m) == BigFloat
    @test m(m.pos) ≈ BigFloat(1)
end

@testset "Gaussian" begin
    test_model_interface(Gaussian)

    @testset "isotropic" begin
        m = Gaussian(fwhm=10)
        expected = exp(-4 * log(2) * sum(abs2, SA[1, 2]) / 100)
        @test m[2, 1] ≈ m(1, 2) ≈ expected
        @test repr(m) == "Gaussian{Float64}(pos=[0, 0], fwhm=10, amp=1.0)"
    end
    
    @testset "diagonal" begin
        m = Gaussian(fwhm=(10, 9))
        wdist = (1/10)^2 + (2/9)^2
        expected = exp(-4 * log(2) * wdist)
        @test m[2, 1] ≈ m(1, 2) ≈ expected
        @test repr(m) == "Gaussian{Float64}(pos=[0, 0], fwhm=(10, 9), amp=1.0)"
    end

    # test Normal alias
    @test Normal(fwhm=10) === Gaussian(fwhm=10)

    @testset "gradients" begin
        FiniteDifferences.to_vec(x::Integer) = Bool[], _ -> x
        # have to make sure PSFs are all floating point so tangents don't have type issues
        psf_iso = Gaussian(fwhm=10.0, pos=zeros(2))
        psf_tang = Tangent{Gaussian}(fwhm=rand(rng), pos=rand(rng, 2), amp=rand(rng), indices=NoTangent())
        point = Float64[1, 2]
        test_frule(psf_iso ⊢ psf_tang, point)
        test_rrule(psf_iso ⊢ psf_tang, point)

        psf_diag = Gaussian(fwhm=Float64[10, 8], pos=zeros(2))
        psf_tang = Tangent{Gaussian}(fwhm=rand(rng, 2), pos=rand(rng, 2), amp=rand(rng), indices=NoTangent())
        test_frule(psf_diag ⊢ psf_tang, point)
        test_rrule(psf_diag ⊢ psf_tang, point)
    end
end


@testset "AiryDisk" begin
    test_model_interface(AiryDisk)

    @testset "isotropic" begin
        m = AiryDisk(fwhm=10)
        radius = m.fwhm * 1.18677
        # first radius is 0
        @test m(radius, 0) ≈ 0 atol=eps(Float64)
        @test m(-radius, 0) ≈ 0 atol=eps(Float64)
        @test m(0, radius) ≈ 0 atol=eps(Float64)
        @test m(0, -radius) ≈ 0 atol=eps(Float64)
        @test repr(m) == "AiryDisk{Float64}(pos=[0, 0], fwhm=10, amp=1.0)"
    end

    @testset "diagonal" begin
        m = AiryDisk(fwhm=(10, 9))
        r1 = m.fwhm[1] * 1.18677
        r2 = m.fwhm[2] * 1.18677
        # first radius is 0
        @test m(r1, 0) ≈ 0 atol=eps(Float64)
        @test m(-r1, 0) ≈ 0 atol=eps(Float64)
        @test m(0, r2) ≈ 0 atol=eps(Float64)
        @test m(0, -r2) ≈ 0 atol=eps(Float64)
        @test repr(m) == "AiryDisk{Float64}(pos=[0, 0], fwhm=(10, 9), amp=1.0)"
    end
end

@testset "Moffat" begin
    test_model_interface(Moffat)

    @testset "isotropic" begin
        m = Moffat(fwhm=10)
        expected = inv(1 + sum(abs2, SA[1, 2]) / 25)
        @test m[2, 1] ≈ m(1, 2) ≈ expected
        @test repr(m) == "Moffat{Float64}(pos=[0, 0], fwhm=10, amp=1.0, alpha=1)"
    end

    @testset "diagonal" begin
        m = Moffat(fwhm=(10, 9))
        wdist = (1/5)^2 + (2/4.5)^2
        expected = inv(1 + wdist)
        @test m[2, 1] ≈ m(1, 2) ≈ expected
        @test repr(m) == "Moffat{Float64}(pos=[0, 0], fwhm=(10, 9), amp=1.0, alpha=1)"
    end

    @testset "alpha" begin
        m = Moffat(fwhm=10, alpha=2)
        expected = inv(1 + sum(abs2, SA[1, 2]) / 25)^2
        @test m[2, 1] ≈ m(1, 2) ≈ expected
        @test repr(m) == "Moffat{Float64}(pos=[0, 0], fwhm=10, amp=1.0, alpha=2)"
    end
end

include("plotting.jl")

using PSFModels: CircularGaussianPSF, GaussianPSF, CircularGaussianPRF, GaussianPRF, evaluate, centroid, integral, evaluate_fg, evaluate_fgh, AbstractPSFModel, extent, render, render!, theta, amplitude, background, fwhm, peak, effective_area
using Test

# Tests generic API, type return, etc
function test_common(model::AbstractPSFModel{T}) where T
    # Evaluation
    @test @inferred(evaluate(model, centroid(model)...)) isa T
    @test @inferred(evaluate(model, CartesianIndex(round.(Int, centroid(model))))) isa T
    @test model(centroid(model)...) ≈ evaluate(model, centroid(model)...)
    @test model(CartesianIndex(round.(Int, centroid(model)))) ≈ evaluate(model, CartesianIndex(round.(Int, centroid(model))))
    ex = @inferred extent(model)
    @test ex isa Tuple{Tuple{T, T}, Tuple{T, T}}
    x, y = range(ex[1][1], ex[1][2]; step=one(T)), range(ex[2][1], ex[2][2]; step=one(T))
    m = evaluate.(model, x, y')
    @test m isa Matrix{T}
    @test size(m) == (length(x), length(y))
    ex = (round.(Int, ex[1]), round.(Int, ex[2]))
    x, y = ex[1][1]:ex[1][2], ex[2][1]:ex[2][2]
    m = evaluate.(model, CartesianIndices((x, y)))
    @test m isa Matrix{T}
    @test size(m) == (length(x), length(y))
    @test @inferred(render(model, (x,y))) ≈ m
    m2 = similar(m)
    fill!(m2, 0)
    @test @inferred(render!(m2, model, (x,y))) ≈ m
    @test m2 ≈ evaluate.(model, x, y')
    @test @inferred(render(model)) isa Matrix{T}

    # API functions
    @test @inferred(centroid(model)) isa Tuple{T, T}
    @test @inferred(integral(model)) isa T
    @test @inferred(peak(model)) isa T
    @test @inferred(amplitude(model)) isa T
    @test @inferred(background(model)) isa T
    # @test effective_area(model) isa T # no generic method yet
    # @test fwhm(model) isa T # no generic method yet
    @test @inferred(theta(model)) isa T
end

for model in (CircularGaussianPSF(x=1.3, y=2.4, fwhm=3.0, flux=120.0, bkg=10.0),
              CircularGaussianPSF(x=1.3f0, y=2.4f0, fwhm=3.0f0, flux=120.0f0, bkg=10.0f0),
              GaussianPSF(x=2.5, y=5.0, x_fwhm=3.0, y_fwhm=4.0, theta=35, flux=120.0, bkg=10),
              GaussianPSF(x=2.5f0, y=5.0f0, x_fwhm=3.0f0, y_fwhm=4.0f0, theta=35f0, flux=120.0f0, bkg=10.0f0),
              CircularGaussianPRF(x=1.3, y=2.4, fwhm=3.0, flux=120.0, bkg=10.0),
              CircularGaussianPRF(x=1.3f0, y=2.4f0, fwhm=3.0f0, flux=120.0f0, bkg=10.0f0),
              GaussianPRF(x=2.5, y=5.0, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=120.0, bkg=10),
              GaussianPRF(x=2.5f0, y=5.0f0, x_fwhm=3.0f0, y_fwhm=4.0f0, theta=35.0f0, flux=120.0f0, bkg=10.0f0))
    @testset "API: $(typeof(model))" begin
        test_common(model)
    end
end

# Test specific models; verify return values
@testset "CircularGaussianPSF" begin
    @testset "constructor promotion" begin
        @test CircularGaussianPSF(x=1.3, y=2.4, fwhm=3.0, flux=120.0, bkg=10) isa CircularGaussianPSF{Float64}
        @test CircularGaussianPSF(x=1.3f0, y=2.4f0, fwhm=3.0f0, flux=120.0f0, bkg=10.0f0) isa CircularGaussianPSF{Float32}
        @test CircularGaussianPSF(x=1, y=2, fwhm=3, flux=120, bkg=10) isa CircularGaussianPSF{Float64}
        @test CircularGaussianPSF(x=BigFloat(1.3), y=BigFloat(2.4), fwhm=BigFloat(3.0), flux=BigFloat(120.0), bkg=BigFloat(10.0)) isa CircularGaussianPSF{BigFloat}
    end

    m = CircularGaussianPSF(x=0, y=0, fwhm=10, flux=1, bkg=10)
    @test centroid(m) == (0.0, 0.0)
    @test integral(m) == 1.0
    @test fwhm(m) == (10.0, 10.0)
    @test effective_area(m) ≈ 226.6180070913597 rtol=1e-6
    @test background(m) == 10.0
    @test peak(m) ≈ amplitude(m) + background(m)
    @test theta(m) == 0.0
    r1 = evaluate(m, 1, 2)
    @test r1 isa Float64
    @test r1 ≈ 10.0076829778398427705 ≈ m(1, 2)
    let (f, g) = evaluate_fg(m, 1, 2)
        @test f ≈ evaluate(m, 1, 2)
        @test g ≈ [0.0004260347542393244, 0.0008520695084786488, -0.001323578190848892, 0.0076829778398427705, 1.0]
    end
    let (f, g, h) = evaluate_fgh(m, 1, 2)
        @test f ≈ evaluate(m, 1, 2)
        @test g ≈ [0.0004260347542393244, 0.0008520695084786488, -0.001323578190848892, 0.0076829778398427705, 1.0]
        @test h ≈ [-0.0004024103711416015 4.724876619544591e-5 -0.00015860171014686828 0.0004260347542393244 0.0; 4.724876619544591e-5 -0.00033153722184843266 -0.00031720342029373656 0.0008520695084786488 0.0; -0.00015860171014686828 -0.00031720342029373656 0.00031777260218123345 -0.0013235781908488918 0.0; 0.0004260347542393244 0.0008520695084786488 -0.0013235781908488918 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
    end
end

@testset "GaussianPSF" begin
    @testset "constructor promotion" begin
        @test GaussianPSF(x=1.3, y=2.4, x_fwhm=3.0, y_fwhm=4.0, theta=35, flux=120.0, bkg=10.0) isa GaussianPSF{Float64}
        @test GaussianPSF(x=1.3f0, y=2.4f0, x_fwhm=3.0f0, y_fwhm=4.0f0, theta=35f0, flux=120.0f0, bkg=10.0f0) isa GaussianPSF{Float32}
        @test GaussianPSF(x=1, y=2, x_fwhm=3, y_fwhm=4, theta=35, flux=120, bkg=10) isa GaussianPSF{Float64}
        @test GaussianPSF(x=BigFloat(1.3), y=BigFloat(2.4), x_fwhm=BigFloat(3.0), y_fwhm=BigFloat(4.0), theta=BigFloat(35), flux=BigFloat(120.0), bkg=BigFloat(10.0)) isa GaussianPSF{BigFloat}
    end

    m = GaussianPSF(x=0, y=0, x_fwhm=10, y_fwhm=6, theta=30, flux=1, bkg=10)
    @test centroid(m) == (0.0, 0.0)
    @test integral(m) == 1.0
    @test fwhm(m) == (10.0, 6.0)
    @test effective_area(m) ≈ 135.9708042548158 rtol=1e-6
    @test background(m) == 10.0
    @test peak(m) ≈ amplitude(m) + background(m)
    @test theta(m) == 30.0
    r1 = evaluate(m, 1, 2)
    @test r1 isa Float64
    @test r1 ≈ 10.011881854589938992
    let (f, g) = evaluate_fg(m, 1, 2)
        @test f ≈ evaluate(m, 1, 2)
        @test g ≈ [-6.269560630626032e-5, 0.0025675279951917654, -0.0009587636050457796, -0.0015172854575571534, 4.7000306666382166e-5, 0.01188185458993899, 1.0]
    end
    let (f, g, h) = evaluate_fgh(m, 1, 2)
        @test f ≈ evaluate(m, 1, 2)
        @test g ≈ [-6.269560630626032e-5, 0.0025675279951917654, -0.0009587636050457796, -0.0015172854575571534, 4.7000306666382166e-5, 0.01188185458993899, 1.0]
        @test h ≈ [-0.0009513701779307769 0.0004936505237653131 -0.00020789110903162168 0.0003838214631391412 -2.986906138973509e-6 -6.26956063062603e-5 0.0; 0.0004936505237653131 -0.0009825507698117963 -0.0003301242561598406 -0.0009787987399725283 -3.547464474559758e-5 0.0025675279951917654 0.0; -0.00020789110903162168 -0.0003301242561598406 0.0001273559772475686 0.00012243190355172497 1.4950134657622455e-6 -0.0009587636050457795 0.0; 0.0003838214631391412 -0.0009787987399725283 0.00012243190355172497 0.0002922935533913823 -3.048115727009258e-5 -0.001517285457557153 0.0; -2.986906138973509e-6 -3.547464474559758e-5 1.4950134657622455e-6 -3.048115727009258e-5 -5.148866586397458e-7 4.700030666638217e-5 0.0; -6.26956063062603e-5 0.0025675279951917654 -0.0009587636050457795 -0.001517285457557153 4.700030666638217e-5 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
    end
    # equal fwhm + theta=0 reduces to CircularGaussianPSF
    mc = CircularGaussianPSF(x=1.5, y=2.5, fwhm=8, flux=3, bkg=0)
    mg = GaussianPSF(x=1.5, y=2.5, x_fwhm=8, y_fwhm=8, theta=0, flux=3, bkg=0)
    @test evaluate(mc, 3, 4) ≈ evaluate(mg, 3, 4)
end

@testset "CircularGaussianPRF" begin
    @testset "constructor promotion" begin
        @test CircularGaussianPRF(x=1.3, y=2.4, fwhm=3.0, flux=120.0, bkg=10) isa CircularGaussianPRF{Float64}
        @test CircularGaussianPRF(x=1.3f0, y=2.4f0, fwhm=3.0f0, flux=120.0f0, bkg=10.0f0) isa CircularGaussianPRF{Float32}
        @test CircularGaussianPRF(x=1, y=2, fwhm=3, flux=120, bkg=10) isa CircularGaussianPRF{Float64}
        @test CircularGaussianPRF(x=BigFloat(1.3), y=BigFloat(2.4), fwhm=BigFloat(3.0), flux=BigFloat(120.0), bkg=BigFloat(10.0)) isa CircularGaussianPRF{BigFloat}
    end

    m = CircularGaussianPRF(x=0, y=0, fwhm=10, flux=1, bkg=10)
    @test centroid(m) == (0.0, 0.0)
    @test integral(m) == 1.0
    @test fwhm(m) == (10.0, 10.0)
    @test effective_area(m) ≈ 226.6180070913597 rtol=1e-6
    @test background(m) == 10.0
    @test peak(m) ≈ amplitude(m) + background(m)
    @test theta(m) == 0.0
    r1 = evaluate(m, 1, 2)
    @test r1 isa Float64
    @test r1 ≈ 10.007652480658708 ≈ m(1, 2)
    let (f, g) = evaluate_fg(m, 1, 2)
        @test f ≈ evaluate(m, 1, 2)
        @test g ≈ [0.00042238646952845753, 0.0008447735386125179, -0.0013132200987741396, 0.007652480658708134, 1.0]
    end
end

@testset "GaussianPRF" begin
    @testset "constructor promotion" begin
        @test GaussianPRF(x=1.3, y=2.4, x_fwhm=3.0, y_fwhm=4.0, theta=35.0, flux=120.0, bkg=10.0) isa GaussianPRF{Float64}
        @test GaussianPRF(x=1.3f0, y=2.4f0, x_fwhm=3.0f0, y_fwhm=4.0f0, theta=35.0f0, flux=120.0f0, bkg=10.0f0) isa GaussianPRF{Float32}
        @test GaussianPRF(x=1, y=2, x_fwhm=3, y_fwhm=4, theta=35, flux=120, bkg=10) isa GaussianPRF{Float64}
        @test GaussianPRF(x=BigFloat(1.3), y=BigFloat(2.4), x_fwhm=BigFloat(3.0), y_fwhm=BigFloat(4.0), theta=BigFloat(35.0), flux=BigFloat(120.0), bkg=BigFloat(10.0)) isa GaussianPRF{BigFloat}
    end

    m = GaussianPRF(x=0, y=0, x_fwhm=10, y_fwhm=6, theta=55.0, flux=1, bkg=10)
    @test centroid(m) == (0.0, 0.0)
    @test integral(m) == 1.0
    @test fwhm(m) == (10.0, 6.0)
    @test effective_area(m) ≈ 135.9708042548158 rtol=1e-6
    @test background(m) == 10.0
    @test peak(m) ≈ amplitude(m) + background(m)
    @test theta(m) == 55.0
    r1 = evaluate(m, 1, 2)
    @test r1 isa Float64
    @test r1 ≈ 10.012636019260277
    let (f, g) = evaluate_fg(m, 1, 2)
        @test f ≈ evaluate(m, 1, 2)
        @test g ≈ [0.00036857849874118024, 0.001625201537357658, -0.0009181257637637815, -0.0020450982870187633, 1.5499301125181995e-5, 0.012636019260278195, 1.0]
    end

    # equal x/y fwhm with theta=0 collapses to CircularGaussianPRF
    mc = CircularGaussianPRF(x=1.5, y=2.5, fwhm=8, flux=3, bkg=0)
    mg = GaussianPRF(x=1.5, y=2.5, x_fwhm=8, y_fwhm=8, theta=0, flux=3, bkg=0)
    @test evaluate(mc, 3, 4) ≈ evaluate(mg, 3, 4)
end
using PSFModels: fit


function generate_model(rng, model, params, inds)
    cartinds = CartesianIndices(inds)
    psf = model.(cartinds; params...)
    return psf
end

function test_fitting(rng, model, params, inds; kwargs...)
    psf = generate_model(rng, model, params, inds)
    _keys = keys(params)
    _vals = PSFModels.vector_from_params(Float64, params)
    # perturb starting guess by a little
    _vals .*= 1 .+ 1e-2 .* randn(rng, length(_vals))
    P0 = PSFModels.generate_params(_keys, _vals)
    P, bestfit = fit(model, P0, psf; kwargs...)
    for k in _keys
        if P[k] isa Tuple
            @test P[k][1] ≈ params[k][1] rtol=1e-2
            @test P[k][2] ≈ params[k][2] rtol=1e-2
        else
            @test P[k] ≈ params[k] rtol=1e-2
        end
    end
    @test bestfit ≈ psf rtol=1e-2
end

@testset "test fitting synthetic PSFs" begin
    inds = 1:30, 1:30
    @testset "gaussian" begin
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds)
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds)
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds)
        
        # test using L1 loss
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds; loss=abs)
        # test using frozen variables
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds; func_kwargs=(;theta=0))

        @test_throws ArgumentError fit(gaussian, (x=12, y=13, fwhm=1, theta=13), randn(10, 10))
    end
    @testset "airydisk" begin
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds)
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=2.6, amp=5, ratio=0.12), inds)
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds)
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds)
        
        # test using L1 loss
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds; loss=abs)
        # test using frozen variables
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds; func_kwargs=(;ratio=0))

        @test_throws ArgumentError fit(airydisk, (x=12, y=13, fwhm=1, theta=13), randn(10, 10))
    end

    @testset "moffat" begin
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=2.6, amp=5, alpha=1.2), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5, alpha=1.2), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5, alpha=1.2), inds)

        # test using L1 loss
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds; loss=abs)
        # test using frozen variables
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds; func_kwargs=(;alpha=1))

        @test_throws ArgumentError fit(moffat, (x=12, y=13, fwhm=1, theta=13), randn(10, 10))
    end
end
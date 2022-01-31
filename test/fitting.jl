using PSFModels: fit


function generate_model(rng, model, params, inds)
    cartinds = CartesianIndices(inds)
    psf = model.(cartinds; params...)
    return psf
end

function test_fitting(rng, model, params, inds; kwargs...)
    psf = generate_model(rng, model, params, inds)
    _keys = keys(params)
    if :fwhm in _keys && params.fwhm isa Tuple
        _ind = findfirst(==(:fwhm), _keys)
        vals = values(params)
        first_half = vals[begin:_ind - 1]
        fwhmx, fwhmy = vals[_ind]
        second_half = vals[_ind + 1:end]
        _vals = Float64[first_half..., fwhmx, fwhmy, second_half...]
    else
        _vals = collect(Float64, values(params))
    end
    # perturb starting guess by a little
    _vals .*= 1 .+ 1e-2 .* randn(rng, length(_vals))
    P, bestfit = fit(model, _keys, _vals, psf; kwargs...)
    for k in _keys
        if P[k] isa Tuple
            @test P[k][1] ≈ params[k][1] rtol=1e-6
            @test P[k][2] ≈ params[k][2] rtol=1e-6
        else
            @test P[k] ≈ params[k] rtol=1e-6
        end
    end
    @test bestfit ≈ psf rtol=1e-6
end

@testset "test fitting synthetic PSFs" begin
    inds = 1:30, 1:30
    @testset "gaussian" begin
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds)
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds)
        test_fitting(rng, gaussian, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds)
    end
    @testset "airydisk" begin
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds)
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=2.6, amp=5, ratio=0.12), inds)
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds)
        test_fitting(rng, airydisk, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds)
    end

    @testset "moffat" begin
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=2.6, amp=5), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=2.6, amp=5, alpha=1.2), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), amp=5, alpha=1.2), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5), inds)
        test_fitting(rng, moffat, (x=13.5, y=12.3, fwhm=(2.6, 2.4), theta=12, amp=5, alpha=1.2), inds)
    end
end
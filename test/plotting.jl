using RecipesBase: apply_recipe

@testset "plotting - $K" for K in (Gaussian, AiryDisk, Moffat)
    psf = K(y=1, fwhm=5)
    recipes = apply_recipe(Dict{Symbol,Any}(), psf)
    for rec in recipes
        @test getfield(rec, 1) == Dict{Symbol,Any}(
            :seriestype => :heatmap,
            :xlim => (-15, 15),
            :ylim => (-14, 16),
            :aspect_ratio => 1

        xs = rec.args[1]
        ys = rec.args[2]
        _psf = rec.args[3]
        
        @test xs == -15:15
        @test ys == -14:16
        @test _psf == collect(psf)
    end
end
using RecipesBase: apply_recipe

@testset "plotting - $K" for K in (Gaussian, AiryDisk, Moffat)
    psf = K(x=0, y=1, fwhm=5)
    recipes = apply_recipe(Dict{Symbol,Any}(), psf)
    for rec in recipes
        @test rec.args[1] === psf
        xs, ys = axes(psf)
        @test rec.args[2] == xs
        @test rec.args[3] == ys
    end

    recipes_full = apply_recipe(Dict{Symbol,Any}(), psf, axes(psf)...)
    for rec in recipes_full
        @test rec.plotattributes == Dict{Symbol,Any}(
            :seriestype => :heatmap,
            :xlims => (-8, 8),
            :ylims => (-7, 9),
            :aspect_ratio => 1,
            :xguide => "x",
            :yguide => "y",
        )

        xs = rec.args[1]
        ys = rec.args[2]
        _psf = rec.args[3]
        
        @test xs == -8:8
        @test ys == -7:9
        @test _psf == transpose(collect(psf))
    end
end
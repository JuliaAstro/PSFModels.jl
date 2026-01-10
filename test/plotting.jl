using PSFModels: PsfPlot
using RecipesBase: apply_recipe

@testset "plotting - $K" for K in (gaussian, airydisk, moffat)
    psf = K(x=0, y=1, fwhm=5)
    inds = (-8:8, -7:9)
    pplot = PsfPlot((psf, inds...))
    recipes = apply_recipe(Dict{Symbol,Any}(), pplot)
    for rec in recipes
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
        @test _psf ≈ transpose(map(psf, CartesianIndices(inds)))
    end

    pplot = PsfPlot((psf, inds))
    recipes_full = apply_recipe(Dict{Symbol,Any}(), pplot)
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
        @test _psf ≈ transpose(map(psf, CartesianIndices(inds)))
    end
end

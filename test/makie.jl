using Makie

const PSFModelsMakieExt = Base.get_extension(PSFModels, :PSFModelsMakieExt)
@test !isnothing(PSFModelsMakieExt)

@testset "makie - $K" for K in (gaussian, airydisk, moffat)
    psf = K(x = 0, y = 1, fwhm = 5)
    xs, ys = -8:8, -7:9
    expected = [float(psf(x, y)) for x in xs, y in ys]

    fig, ax, plt = psfplot(psf, xs, ys)
    @test plt isa PSFModelsMakieExt.PsfPlot
    @test plt.psfdata[] ≈ expected
    # axis defaults matching the Plots.jl recipe
    @test ax.aspect[] == Makie.DataAspect()
    @test ax.xlabel[] == "x"
    @test ax.ylabel[] == "y"
    # the recipe supports colorbars
    Makie.Colorbar(fig[1, 2], plt)

    # tuple form unpacks the indices
    _, _, plt2 = psfplot(psf, (xs, ys))
    @test plt2.psfdata[] ≈ expected

    # mutating form
    fig3 = Makie.Figure()
    ax3 = Makie.Axis(fig3[1, 1])
    plt3 = psfplot!(ax3, psf, xs, ys)
    @test plt3.psfdata[] ≈ expected
end

@testset "makie psfplotview - $K" for K in (gaussian, airydisk, moffat)
    psf = K(x = 0, y = 1, fwhm = 5)
    xs, ys = -8:8, -7:9
    expected = [float(psf(x, y)) for x in xs, y in ys]

    fig, view = psfplotview(psf, xs, ys)
    @test view isa PSFModelsMakieExt.PsfPlotView
    @test view.plt.psfdata[] ≈ expected
    @test view.ax.aspect[] == Makie.DataAspect()
    @test view.ax.xlabel[] == "x"
    @test view.ax.ylabel[] == "y"
    # axis + colorbar
    @test count(x -> x isa Makie.Axis, view.blocks) == 1
    @test count(x -> x isa Makie.Colorbar, view.blocks) == 1

    # tuple form, axis overrides, no colorbar, placed in an existing figure
    fig2 = Makie.Figure()
    view2 = psfplotview(
        fig2[1, 1], psf, (xs, ys);
        colorbar = false, axis = (; title = "PSF", xlabel = "u"),
    )
    @test view2.plt.psfdata[] ≈ expected
    @test view2.ax.title[] == "PSF"
    @test view2.ax.xlabel[] == "u"
    @test count(x -> x isa Makie.Colorbar, view2.blocks) == 0

    # colormapping attributes pass through to the plot
    _, view3 = psfplotview(psf, xs, ys; colorscale = log10, colormap = :inferno)
    @test view3.plt.colorscale[] === log10
end

@testset "makie psfplotview layout" begin
    psf = gaussian(x = 0, y = 1, fwhm = 5)
    xs, ys = -8:8, -7:9   # 17x17 cells -> data aspect ratio 1

    # standalone view: figure is shrink-wrapped onto the aspect-locked panel
    fig, view = psfplotview(psf, xs, ys)
    @test Tuple(Makie.widths(Makie.viewport(fig.scene)[])) != (600, 450)
    layout = view.layout
    @test (layout.colsizes[1] isa Makie.GridLayoutBase.Aspect) ||
        (layout.rowsizes[1] isa Makie.GridLayoutBase.Aspect)

    # the axis box has the data aspect, and the colorbar sits flush against it
    Makie.update_state_before_display!(fig)
    axbox = view.ax.layoutobservables.computedbbox[]
    w, h = Makie.widths(axbox)
    @test w / h ≈ 1.0 atol = 0.05
    cb = only(filter(x -> x isa Makie.Colorbar, view.blocks))
    gap = Makie.left(cb.layoutobservables.computedbbox[]) - Makie.right(axbox)
    @test 0 <= gap <= 40

    # wide data extents lock to the right ratio
    figw, vieww = psfplotview(psf, -30:30, ys)
    Makie.update_state_before_display!(figw)
    ww, hw = Makie.widths(vieww.ax.layoutobservables.computedbbox[])
    @test ww / hw ≈ 61 / 17 rtol = 0.05

    # an explicit figure size wins over the shrink-wrap
    fig2, _ = psfplotview(psf, xs, ys; figure = (; size = (321, 322)))
    @test Tuple(Makie.widths(Makie.viewport(fig2.scene)[])) == (321, 322)

    # an explicitly sized axis derives the other dimension from the data
    # aspect and shrink-wraps via resize_to_layout!
    fig3, view3 = psfplotview(psf, xs, ys; axis = (; height = 200))
    @test view3.ax.width[] ≈ 200.0
    @test all(!isnothing, view3.layoutobservables.autosize[])
    figsize3 = Makie.widths(Makie.viewport(fig3.scene)[])
    @test 200 < figsize3[2] < 320   # 200 + decorations, not the 450 default

    # placed into an existing figure, the parent is never resized
    fig4 = Makie.Figure(size = (500, 500))
    psfplotview(fig4[1, 1], psf, xs, ys)
    @test Tuple(Makie.widths(Makie.viewport(fig4.scene)[])) == (500, 500)
end

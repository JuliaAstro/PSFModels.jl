using PSFModels: GaussianPSF, CircularGaussianPSF, evaluate, _make_fgh, free_params, model_from_vector, fit, fit_lm, render, TukeyLoss, bicubic_interpolate, _fill_missing_bicubic!
using BenchmarkTools
import LossFunctions
import Optim
using PrettyTables: pretty_table

function show_benchmarks(results)
    # Collect results
    sorted  = sort(collect(results), by=first)
    names   = [k for (k,_) in sorted]
    trials  = [v for (_,v) in sorted]

    # Pack into matrix
    data = hcat(
        names,
        [BenchmarkTools.prettytime(median(t).time) for t in trials],
        [BenchmarkTools.prettymemory(median(t).memory) for t in trials],
        [median(t).allocs for t in trials]
    )

    # Make pretty table
    pretty_table(data;
        column_labels = ["Benchmark", "Median Time", "Memory", "Allocs"],
        alignment     = [:l, :r, :r, :r]
    )
end

const SUITE = BenchmarkGroup()
SUITE["parametric"] = BenchmarkGroup()

let model = CircularGaussianPSF(x=15.0, y=15.0, fwhm=4.0, flux=10.0, bkg=1.0)
    inds  = (1:30, 1:30)
    image = render(model, inds)
    fixed = (; bkg=1.1,)
    free_names, free_idx, x0 = free_params(model, fixed)
    free_names = Val(free_names)
    SUITE["parametric"]["circular_gaussian_objective"] = @benchmarkable fgh.fgh(true, G, H, $x0) setup = (
        fgh = _make_fgh($model, $free_names, $free_idx, $fixed, $image, $inds, $nothing, $(LossFunctions.L2DistLoss()));
                        G = zeros(5); H = zeros(5, 5))
end

# ---------------------------------------------------------------------------
# LM fitting benchmarks
# ---------------------------------------------------------------------------
SUITE["fitting"] = BenchmarkGroup()

let model = CircularGaussianPSF(x=15.0, y=15.0, fwhm=4.0, flux=10.0, bkg=1.0)
    inds  = (1:30, 1:30)
    image = render(model, inds)
    init = CircularGaussianPSF(x=15.5, y=14.5, fwhm=3.5, flux=9.0, bkg=1.2)
    SUITE["fitting"]["fit_lm (L2)"] = @benchmarkable fit_lm($init, $image, $inds)
    SUITE["fitting"]["fit_lm (Huber IRLS)"] = @benchmarkable fit_lm($init, $image, $inds;
        reweight=$(LossFunctions.HuberLoss(1.0)))
    SUITE["fitting"]["fit_lm (Tukey IRLS)"] = @benchmarkable fit_lm($init, $image, $inds;
        reweight=$(TukeyLoss()))
    SUITE["fitting"]["fit (Optim NewtonTrustRegion)"] = @benchmarkable fit($init, $image, $inds;
        x_abstol=1e-8, alg=Optim.NewtonTrustRegion())
    SUITE["fitting"]["fit (Optim LBFGS)"] = @benchmarkable fit($init, $image, $inds;
        x_abstol=1e-8, alg=Optim.LBFGS())
end

# ---------------------------------------------------------------------------
# Empirical model benchmarks
# ---------------------------------------------------------------------------
SUITE["empirical"] = BenchmarkGroup()
SUITE["empirical"]["bicubic_interpolate"] = @benchmarkable bicubic_interpolate(x, $3.5, $3.5) setup=(x=rand(7,7))
SUITE["empirical"]["_fill_missing_bicubic!"] = @benchmarkable _fill_missing_bicubic!(x) setup=(x=rand(21,21); inds=([9, 4, 15, 13, 1, 1], [6, 15, 17, 1, 17, 1]); x[inds...] .= NaN) evals=1

for n in (50, 100)
    SUITE["empirical"]["ImagePSF fit, n=$n"] = @benchmarkable psf, result =
        PSFModels.fit(ImagePSF, image, sources.x, sources.y;
            fit_rad = 5.0, oversampling = 2, smooth = true, recenter = false,
            reweight = nothing) setup=(begin
                truth_model = CircularGaussianPRF(x = 0, y = 0, fwhm = 1.8, flux = 1, bkg = 0)
                image, sources = simulate_image((96, 96), truth_model, $n;
                    background = 20.0, noise = :none, flux = (600.0, 900.0),
                    min_separation = 7, border = 8, model_radius = 6)
            end) evals=1
end

# If not on CI, we'll show a nice table
if get(ENV, "CI", "false") == "false"
    # Run the benchmarks
    results = run(SUITE, verbose=true)
    println("⎯⎯⎯ Parametric Suite ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯")
    show_benchmarks(results["parametric"])
    println("⎯⎯⎯ Fitting Suite ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯")
    show_benchmarks(results["fitting"])
    println("⎯⎯⎯ Empirical Suite ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯")
    show_benchmarks(results["empirical"])
end

using PSFModels: GaussianPSF, CircularGaussianPSF, evaluate, _make_fgh, free_params, model_from_vector, fit
using BenchmarkTools
import LossFunctions
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
SUITE["core"] = BenchmarkGroup()

let model = CircularGaussianPSF(x=15.0, y=15.0, fwhm=4.0, flux=10.0, bkg=1.0),
    inds  = (1:30, 1:30),
    image = [evaluate(model, px, py) for px in inds[1], py in inds[2]]
    fixed = (; bkg=1.1,)
    free_names, free_idx, x0 = free_params(model, fixed)
    free_names = Val(free_names)
    SUITE["core"]["circular_gaussian_objective"] = @benchmarkable fgh.fgh(true, G, H, $x0) setup = (
        fgh = _make_fgh($model, $free_names, $free_idx, $fixed, $image, $inds, $nothing, $(LossFunctions.L2DistLoss()));
                        G = zeros(5); H = zeros(5, 5))
end

# If not on CI, we'll show a nice table
if get(ENV, "CI", "false") == "false"
    # Run the benchmarks
    results = run(SUITE, verbose=true)
    println("⎯⎯⎯ Core Suite ⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯⎯")
    show_benchmarks(results["core"])
end

# Benchmarks

The benchmarks can be found in the [`bench/`](https://github.com/JuliaAstro/PSFModels.jl/tree/main/bench) folder. To run them, first install the python dependencies

```
$ cd bench
$ poetry install
$ poetry shell
```
then get the Julia project set up
```
$ PYTHON=$(which python) julia --project=@. -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
```

Then run the benchmark

```
$ julia --project=. bench.jl
```

**System Information**

```
Julia Version 1.6.0
Commit f9720dc2eb* (2021-03-24 12:55 UTC)
Platform Info:
  OS: macOS (x86_64-apple-darwin20.3.0)
  CPU: Intel(R) Core(TM) i5-8259U CPU @ 2.30GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-11.0.1 (ORCJIT, skylake)
Environment:
  JULIA_NUM_THREADS = 1
```

---

### Evaluation benchmark

This benchmark tests how long it takes to evaluate a single point in the PSF model. This may seem contrived, but we expect performance to scale directly from this measure: if it takes 1 microsecond to evaluate a single point, it should take ~1 second to evaluate a 1000x1000 image, with speedups potentially from multithreading or SIMD loop evaluation.

```@setup bench
using CSV, DataFrames
using StatsPlots
benchdir(args...) = joinpath("..", ".." ,"bench", args...);
```


```@example bench
table = CSV.File(benchdir("results.csv")) |> DataFrame
```

```@example bench
@df table groupedbar(:name, [:psfmodels :astropy];
    ylabel="time (s)", yscale=:log10, leg=:outertopright,
    label=["PSFModels.jl" "Astropy"], size=(500, 300))
```

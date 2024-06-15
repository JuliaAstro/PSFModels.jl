# Benchmarks

The benchmarks can be found in the [`bench/`](https://github.com/JuliaAstro/PSFModels.jl/tree/main/bench) folder. To run them, first instantiate the environment

```sh
$ julia --project=bench -e "using Pkg; Pkg.instantiate()"
```

then execute the `bench/bench.jl` file

```sh
$ julia --project=bench bench/bench.jl
```

**System Information**

```plain
Julia Version 1.8.0-DEV.1437
Commit a0093d2ffb (2022-02-01 00:11 UTC)
Platform Info:
  OS: macOS (arm64-apple-darwin21.2.0)
  CPU: Apple M1 Max
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-13.0.0 (ORCJIT, cyclone)
Environment:
  JULIA_NUM_THREADS = 1
```

---

### Evaluation benchmark

This benchmark tests how long it takes to evaluate a single point in the PSF model. This may seem contrived, but we expect performance to scale directly from this measure: if it takes 1 microsecond to evaluate a single point, it should take ~1 second to evaluate a 1000×1000 image, with speedups potentially from multithreading or SIMD loop evaluation.

```@setup bench
using CSV
using DataFrames
using StatsPlots
benchdir(args...) = joinpath("..", ".." ,"bench", args...);
```


```@example bench
table = CSV.read(benchdir("evaluation_results.csv"), DataFrame)
show(table) # hide
```

```@example bench
@df table groupedbar(
    :name, [:psfmodels :astropy];
    ylabel="time (s)", yscale=:log10, legend=:outertopright,
    label=["PSFModels.jl" "Astropy"], size=(500, 300),
)
```


### Fitting benchmark

This benchmark tests how long it takes to fit a PSF Model to a stamp with size (39, 39). In all cases, we use equivalently complex models, the default fitters for PSFModels.jl, and the default `LevMarLSQFit` fitter for astropy.

```@example bench
table = CSV.read(benchdir("fitting_results.csv"), DataFrame)
show(table) # hide
```

```@example bench
@df table groupedbar(
    :name, [:psfmodels :astropy];
    ylabel="time (s)", yscale=:log10, legend=:outertopright,
    label=["PSFModels.jl" "Astropy"], size=(500, 300),
)
```

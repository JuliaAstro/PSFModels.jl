```@meta
CurrentModule = PSFModels
```

# PSFModels.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/juliaastro/PSFModels.jl)
[![Build Status](https://github.com/juliaastro/PSFModels.jl/workflows/CI/badge.svg?branch=main)](https://github.com/juliaastro/PSFModels.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/P/PSFModels.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![Coverage](https://codecov.io/gh/juliaastro/PSFModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/juliaastro/PSFModels.jl)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Installation

PSFModels can be added from the Julia package manager

```julia
julia>]

(@v1.6) pkg> add PSFModels
```

## Getting Started

To import the library

```julia
julia> using PSFModels
```

None of the models are exported to avoid namespace clashes, but it can be verbose to continuously rewrite `PSFModels`. You can either import names directly

```julia
julia> using PSFModels: Gaussian

julia> model = Gaussian(fwhm=8)
```

or you can create an alias for `PSFModels`

```julia
# julia version 1.5 or below
using PSFModels
const M = PSFModels
# julia version 1.6 or above
import PSFModels as M

model = M.Gaussian(fwhm=10)
```

```@docs
PSFModels
```

## Benchmarks

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

### Evaluation benchmark

This benchmark tests how long it takes to evaluate a single point in the PSF model. This may seem contrived, but we expect performance to scale directly from this measure: if it takes 1 microsecond to evaluate a single point, it should take ~1 second to evaluate a 1000x1000 image, with speedups potentially from multithreading or SIMD loop evaluation.

```@setup bench
using CSV, DataFrames
using StatsPlots
benchdir(args...) = joinpath("..", ".." ,"bench", args...);
```

---

```@example bench
table = CSV.File(benchdir("results.csv")) |> DataFrame
```

```@example bench
@df table groupedbar(:name, [:psfmodels :astropy];
    ylabel="time (s)", yscale=:log10, leg=:outertopright,
    label=["PSFModels.jl" "Astropy"], size=(500, 300))
```

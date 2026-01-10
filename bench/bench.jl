using BenchmarkTools
using CSV
using DataFrames
using HCIDatasets
using PSFModels
using PythonCall

@info "Starting PSF evaluation benchmark"

ap = pyimport("astropy")
apm = pyimport("astropy.modeling.models")
apf = pyimport("astropy.modeling.fitting")

@info "Using astropy version $(ap.__version__)"


julia_models = [
    gaussian(x=0, y=0, fwhm=(5.3, 4.7)),
    airydisk(x=0, y=0, fwhm=12.7),
    moffat(x=0, y=0, fwhm=8.1)
]
astropy_models = [
    apm.Gaussian2D(x_stddev=5.3, y_stddev=4.7),
    apm.AiryDisk2D(radius=12.7),
    apm.Moffat2D(gamma=8.1, alpha=-1)
]

names = ["Gaussian", "AiryDisk", "Moffat"]
jl_times = []
py_times = []

for (jl_mod, py_mod, name) in zip(julia_models, astropy_models, names)
    jl_time = @belapsed $jl_mod(0.5, 0.7)
    push!(jl_times, jl_time)

    py_time = @belapsed $py_mod(0.5, 0.7)
    push!(py_times, py_time)

    @info "$name" PSFModels=jl_time astropy=py_time
end

filename = joinpath(@__DIR__, "evaluation_results.csv")
DataFrame(name=names, psfmodels=jl_times, astropy=py_times) |> CSV.write(filename)

@info "Results saved to $filename"

@info "Finished PSF evaluation benchmark"

################################################################################
################################################################################
################################################################################

@info "Starting PSF fitting benchmark"

psf = HCIDatasets.BetaPictoris[:psf]
jl_times = []
py_times = []


P0 = (x=20, y=20, fwhm=(5, 5), amp=0.1)
t_jl = @belapsed PSFModels.fit(gaussian, $P0, $psf)
cartinds = CartesianIndices(psf)
axx = map(idx -> idx.I[1], cartinds)
axy = map(idx -> idx.I[2], cartinds)
t_py = @belapsed begin
    g0 = $apm.Gaussian2D(x_stddev=5/2.355, y_stddev=5/2.355, x_mean=20, y_mean=20, amplitude=0.1)
    fitter = $apf.LevMarLSQFitter()
    fitter(g0, $axx, $axy, $psf)
end
push!(jl_times, t_jl)
push!(py_times, t_py)
@info "Gaussian" PSFModels=t_jl astropy=t_py


P0 = (x=20, y=20, fwhm=5, amp=0.1)
t_jl = @belapsed PSFModels.fit(airydisk, $P0, $psf)
t_py = @belapsed begin
    a0 = $apm.AiryDisk2D(radius=5/0.973, x_0=20, y_0=20, amplitude=0.1)
    fitter = $apf.LevMarLSQFitter()
    fitter(a0, $axx, $axy, $psf)
end
push!(jl_times, t_jl)
push!(py_times, t_py)
@info "AiryDisk" PSFModels=t_jl astropy=t_py

P0 = (x=20, y=20, fwhm=5, amp=0.1)
t_jl = @belapsed PSFModels.fit(moffat, $P0, $psf)
t_py = @belapsed begin
    m0 = $apm.Moffat2D(gamma=5, x_0=20, y_0=20, amplitude=0.1)
    fitter = $apf.LevMarLSQFitter()
    fitter(m0, $axx, $axy, $psf)
end
push!(jl_times, t_jl)
push!(py_times, t_py)
@info "Moffat" PSFModels=t_jl astropy=t_py

filename = joinpath(@__DIR__, "fitting_results.csv")
DataFrame(name=names, psfmodels=jl_times, astropy=py_times) |> CSV.write(filename)

@info "Results saved to $filename"

@info "Finished PSF fitting benchmark"

################################################################################
################################################################################
################################################################################

nothing

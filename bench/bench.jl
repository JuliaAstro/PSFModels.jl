using BenchmarkTools
using CSV
using DataFrames
using PSFModels
using PythonCall

@info "Starting PSF evaluation benchmark"

ap = pyimport("astropy")
apm = pyimport("astropy.modeling.models")
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

nothing

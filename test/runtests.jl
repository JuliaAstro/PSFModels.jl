using SafeTestsets

@safetestset "API tests" include("common_tests.jl")
@safetestset "Functional interface - API" include("functional_test.jl")
@safetestset "Functional interface - Fitting" include("functional_fitting_test.jl")
@safetestset "Fitting" include("fitting_struct.jl")
@safetestset "Plotting" include("plotting.jl")

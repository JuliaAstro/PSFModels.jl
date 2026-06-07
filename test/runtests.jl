using SafeTestsets

@safetestset "API tests" include("common_tests.jl")
@safetestset "Functional interface - API" include("functional_test.jl")
@safetestset "Functional interface - Fitting" include("functional_fitting_test.jl")
@safetestset "Fitting" include("fitting_struct.jl")
@safetestset "Empirical models" include("empirical_model_tests.jl")
@safetestset "Plotting" include("plotting.jl")
@safetestset "Simulation" include("simulation_test.jl")
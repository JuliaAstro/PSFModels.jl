
"""
    ScaledPSFModel(amp, model) <: PSFModel

A lazy wrapper for a [`PSFModel`](@ref) that adds the scalar amplitude `amp` to the given `model`. This is convenient to avoid broadcasting and materializing an array from a PSFModel when the scalar multiplication can occur during the function call.

# Examples
```jldoctest
julia> g = PSFModels.Gaussian(10);

julia> g(0, 0)
1.0

julia> m1 = 20 * g; # construct implicitly using scalar * or /

julia> m1(0, 0)
20.0

julia> maximum(g / 100) # scalar division works, too
0.01

julia> m2 = PSFModels.ScaledPSFModel(20, g); # construct explicitly

julia> m1 ≈ m2
true
```
"""
struct ScaledPSFModel{T,M<:PSFModel} <: PSFModel{T}
    amp::T
    model::M
    function ScaledPSFModel(amp::T, m::PSFModel{S}) where {T,S}
        V = promote_type(T,S)
        new{V,typeof(m)}(convert(V, amp), m)
    end
end

Base.size(m::ScaledPSFModel) = size(m.model)
Base.axes(m::ScaledPSFModel) = axes(m.model)

(m::ScaledPSFModel{T})(point::AbstractVector) where {T} = convert(T, m.amp * m.model(point))

Base.:(*)(a::Number, m::PSFModel) = ScaledPSFModel(a, m)
Base.:(*)(m::PSFModel, a::Number) = *(a, m)
Base.:(*)(a::Number, m::ScaledPSFModel) = ScaledPSFModel(a * m.amp, m.model)
Base.:(*)(m::ScaledPSFModel, a::Number) = *(a, m)

Base.:(/)(m::PSFModel, a::Number) = ScaledPSFModel(1 / a, m)
Base.:(/)(m::ScaledPSFModel, a::Number) = ScaledPSFModel(m.amp / a, m)

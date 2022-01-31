using Optim
using Statistics

Model = Union{typeof(gaussian), typeof(normal), typeof(airydisk), typeof(moffat)}

function fit(model::Model, params, X0, image::AbstractMatrix{T}, inds=axes(image); alg=LBFGS(), kwargs...) where T
    if length(params) == length(X0) && :theta in params
        throw(ArgumentError("cannot fit theta for isotropic PSF"))
   end
    cartinds = CartesianIndices(inds)
    function loss(X::AbstractVector{T}) where T
        P = generate_params(params, X)
        minind = map(minimum, inds)
        maxind = map(maximum, inds)
        minind[1] - 0.5 ≤ P.x ≤ maxind[1] + 0.5 || return T(Inf)
        minind[2] - 0.5 ≤ P.y ≤ maxind[2] + 0.5 || return T(Inf)
        all(>(0), P.fwhm) || return T(Inf)
        if :ratio in params
            0 < P.ratio < 1 || return T(Inf)
        end
        if :theta in params
            -45 < P.theta < 45 || return T(Inf)
        end
        # mean square error
        mse = mean(cartinds) do idx
            resid = model(T, idx; P...) - image[idx]
            return resid^2
        end
        return mse
    end
    result = optimize(loss, X0, alg; autodiff=:forward, kwargs...)
    Optim.converged(result) || @warn "optimizer did not converge" result
    X = Optim.minimizer(result)
    P_best = generate_params(params, X)
    return P_best, model.(T, cartinds; P_best...)
end

function generate_params(names, values)
    if length(values) > length(names) && :fwhm in names
        _ind = findfirst(==(:fwhm), names)
        fwhm = values[_ind], values[_ind + 1]
        first_half = @views zip(names[begin:_ind - 1], values[begin:_ind - 1])
        second_half = @views zip(names[_ind + 1:end], values[_ind + 2:end])
        return (;first_half..., fwhm, second_half...)        
    end
    return (;zip(names, values)...)
end

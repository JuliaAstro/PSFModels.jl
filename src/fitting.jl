
const Model = Union{typeof(gaussian), typeof(normal), typeof(airydisk), typeof(moffat)}

"""
    PSFModels.fit(model, params, X0, image, inds=axes(image); alg=LBFGS(), kwargs...)

Fit a PSF model (`model`) defined by the given `params` as a collection of symbols, and initial guess vector `X0`, to the data in `image` at the specified `inds` (by default, the entire array). Additional keyword arguments, as well as the fitting algorithm `alg`, are passed to `Optim.optimize`. By default we use forward-mode autodifferentiation (AD) to derive Jacobians for the fitting functions. Refer to the [Optim.jl documentation](https://julianlsolvers.github.io/Optim.jl/stable/) for more information.

# Choosing parameters

The `fit` function is very powerful because it gives you a great amount of flexibility in the way you fit your models. To demonstrate this, let's start with a simple isotropic [`gaussian`](@ref).

```julia
model = gaussian
# match params to arguments of PSF
params = (:x, :y, :fwhm, :amp)
# initial-guess vector follows order of `params`
X0 = [20, 20, 3, 1]
```

Note that `params` can follow any order

```julia
params = (:amp, :x, :y, :fwhm)
X0 = [1, 20, 20, 3]
```

Now, to extend this interface to the bivariate PSF case, where `fwhm` is two values, all you need to do is add an additional entry in the appropriate position in the intial guess vector

```julia
params = (:x, :y, :fwhm)
# x, y, fwhmx, fwhmy
X0 = [20, 20, 3, 3]
```

and, again, the order does not matter

```julia
model = moffat
params = (:alpha, :x, :y, :fwhm, :amp)
# alpha, x, y, fwhmx, fwhmy, amp
X0 = [1, 20, 20, 3, 3, 10]
```
    
As of right now, the only thing you cannot do is "freeze" a parameter (to do that will require writing your own fitting methods and loss functions).

# Fitting a PSF

After selecting your model and parameters, fitting data is easy

```julia
P = (x=12, y=13, fwhm=13.2, amp=0.1)
psf = gaussian.(CartesianIndicies(1:25, 1:15); P...)

params, synthpsf = PSFModels.fit(gaussian, keys(P), values(P), psf)
```

here `params` is a named tuple of the best fitting parameters. This allows you to instantly create your best-fitting model like

```julia
model = gaussian(; params...)
```

`synthpsf` is the best-fitting model evaluated on the input grid, for direct comparison with the input data.
"""
function fit(model::Model, params, X0::AbstractVector, image::AbstractMatrix{T}, inds=axes(image); alg=LBFGS(), kwargs...) where T
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

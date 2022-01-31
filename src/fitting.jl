
const Model = Union{typeof(gaussian), typeof(normal), typeof(airydisk), typeof(moffat)}

"""
    PSFModels.fit(model, params, image, inds=axes(image); func_kwargs=(;), alg=LBFGS(), kwargs...)

Fit a PSF model (`model`) defined by the given `params` as a named tuple of the parameters to fit and their default values. This model is fit to the data in `image` at the specified `inds` (by default, the entire array). To pass extra keyword arguments to the `model` (i.e., to "freeze" a parameter), pass them in a named tuple to `func_kwargs`. Additional keyword arguments, as well as the fitting algorithm `alg`, are passed to `Optim.optimize`. By default we use forward-mode auto-differentiation (AD) to derive Jacobians for the LBFGS optimization algorithm. Refer to the [Optim.jl documentation](https://julianlsolvers.github.io/Optim.jl/stable/) for more information.

# Choosing parameters

The `fit` function is very powerful because it gives you a great amount of flexibility in the way you fit your models. To demonstrate this, let's start with a simple isotropic [`gaussian`](@ref).

```julia
model = gaussian
# match params to arguments of PSF
params = (x=20, y=20, fwhm=3, amp=1)
```

Note that `params` can follow any order

```julia
params = (amp=1, x=20, y=20, fwhm=3)
```

Now, to extend this interface to the bivariate PSF case, where `fwhm` is two values, all you need to do is use a tuple or vector

```julia
params = (x=20, y=20, fwhm=(3, 3))
```

and, again, the order does not matter

```julia
model = moffat
params = (alpha=1, x=20, y=20, fwhm=3, amp=10)
```

# Fitting a PSF

After selecting your model and parameters, fitting data is easy

```julia
P = (x=12, y=13, fwhm=13.2, amp=0.1)
psf = gaussian.(CartesianIndicies(1:25, 1:15); P...)

params, synthpsf = PSFModels.fit(gaussian, P, psf)
```

here `params` is a named tuple of the best fitting parameters. This allows you to instantly create your best-fitting model like

```julia
model = gaussian(; params...)
```

`synthpsf` is the best-fitting model evaluated on the input grid, for direct comparison with the input data.

To alter parameters without fitting them (i.e., "freeze" them) use `func_kwargs`

```julia
P = (x=12, y=13, fwhm=(12.4, 13.2), amp=0.1)
func_kwargs = (alpha=2)
params, synthpsf = PSFModels.fit(moffat, P, psf; func_kwargs)
```
"""
function fit(model::Model, 
             params, 
             image::AbstractMatrix{T}, 
             inds=axes(image); 
             func_kwargs=(;), 
             alg=LBFGS(), 
             kwargs...) where T
    _keys = keys(params)
    cartinds = CartesianIndices(inds)
    function loss(X::AbstractVector{T}) where T
        P = generate_params(_keys, X)
        minind = map(minimum, inds)
        maxind = map(maximum, inds)
        minind[1] - 0.5 ≤ P.x ≤ maxind[1] + 0.5 || return T(Inf)
        minind[2] - 0.5 ≤ P.y ≤ maxind[2] + 0.5 || return T(Inf)
        all(>(0), P.fwhm) || return T(Inf)
        if :ratio in _keys
            0 < P.ratio < 1 || return T(Inf)
        end
        if :theta in _keys
            -45 < P.theta < 45 || return T(Inf)
        end
        # mean square error
        mse = mean(cartinds) do idx
            resid = model(T, idx; P..., func_kwargs...) - image[idx]
            return resid^2
        end
        return mse
    end
    X0 = vector_from_params(T, params)
    result = optimize(loss, X0, alg; autodiff=:forward, kwargs...)
    Optim.converged(result) || @warn "optimizer did not converge" result
    X = Optim.minimizer(result)
    P_best = generate_params(_keys, X)
    return P_best, model.(T, cartinds; P_best..., func_kwargs...)
end

function vector_from_params(T, params)
    _keys = keys(params)
    _vals = values(params)
    if :fwhm in _keys && params.fwhm isa BivariateLike
        _ind = findfirst(==(:fwhm), _keys)
        first_half = @views _vals[begin:_ind - 1]
        fwhmx, fwhmy = _vals[_ind]
        second_half = @views _vals[_ind + 1:end]
        return T[first_half..., fwhmx, fwhmy, second_half...]
    end
    if :theta in _keys && !(:fwhm isa BivariateLike)
        throw(ArgumentError("cannot fit theta for isotropic distribution"))
    end
    return collect(T, _vals)
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


const Model = Union{typeof(gaussian), typeof(normal), typeof(airydisk), typeof(moffat)}

"""
    PSFModels.fit(model, params, image, inds=axes(image);
                  func_kwargs=(;), loss=abs2, maxfwhm=Inf, 
                  alg=Optim.NewtonTrustRegion(), inv_var=nothing,
                  kwargs...)

Fit the parameters of a PSF model (`model`) with initial guess `params`.
This model is fit to the data in `image` at the specified `inds` 
(by default, the entire array). Inverse variance
data weights can be passed via the `inv_var` keyword argument, 
which should be the same size as `image` if provided.
The best-fitting parameters are returned as a named tuple, along with the 
best-fitting model. If `inv_var` is provided, the covariance matrix of the 
best-fitting parameters is also returned.

To pass extra keyword arguments to the `model` (i.e., to "freeze" a parameter), 
pass them in a named tuple to `func_kwargs`. The default loss function is the 
chi-squared loss, which uses the the square of the difference (i.e., the L2 norm). 
You can change this to the L1 norm, for example, by passing `loss=abs`. 
The maximum FWHM can be set with `maxfwhm` as a number or tuple.

Additional keyword arguments, as well as the fitting algorithm `alg`, are passed
to `Optim.optimize` as an `Optim.Option`. By default we use forward-mode auto-differentiation (AD) to
derive Jacobians for the
[Newton with Trust Region](https://julianlsolvers.github.io/Optim.jl/stable/#algo/newton_trust_region/)
optimization algorithm. Refer to the
[Optim.jl documentation](https://julianlsolvers.github.io/Optim.jl/stable/)
for more information.

# Arguments
- `model::Model`: the PSF model to fit, e.g., `gaussian` or `moffat`
- `params`: a named tuple of the parameters to fit and their default values, \
e.g., `(x=20, y=20, fwhm=3)`
- `image::AbstractMatrix{T}`: the data to fit the model to
- `inds`: the indices of the data to fit to, by default the entire array

# Keyword arguments
- `func_kwargs`: a named tuple of extra keyword arguments to pass to the model; \
this enables freezing a parameter, e.g., `(fwhm=3)`
- `loss`: the loss function to minimize, by default `abs2`, the chi-squared loss (L2 norm)
- `maxfwhm`: the maximum FWHM, can be a number or tuple
- `alg`: the optimization algorithm to use, by default `Optim.NewtonTrustRegion()`
- `inv_var`: inverse variance weights on the input `image`; if provided, \
must be the same size as `image`, and the covariance matrix of the best-fitting \
parameters will be returned
- `kwargs`: additional keyword arguments passed to `Optim.Options`

# Returns
- `P_best`: a named tuple of the best-fitting parameters
- `best_model`: the best-fitting model, for direct comparison with the input data
If `inv_var` is given, also returns:
- `cov`: the covariance matrix of the best-fitting parameters. Standard errors on the parameters can be estimated as `sqrt.(LinearAlgebra.diag(cov))`.

# Choosing parameters

The `fit` function is very powerful because it gives you a great amount of
flexibility in the way you fit your models. To demonstrate this, let's start
with a simple isotropic [`gaussian`](@ref).

```julia
model = gaussian
# match params to arguments of PSF
params = (x=20, y=20, fwhm=3, amp=1)
```

Note that `params` can follow any order

```julia
params = (amp=1, x=20, y=20, fwhm=3)
```

Now, to extend this interface to the bivariate PSF case, where `fwhm` is two
values, all you need to do is use a tuple or vector

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
psf = gaussian.(CartesianIndices((1:25, 1:15)); P...)

params, synthpsf = PSFModels.fit(gaussian, P, psf)
```

here `params` is a named tuple of the best fitting parameters.
It will not include any fixed parameters.

`synthpsf` is the best-fitting model, for direct comparison with the input data.

```julia
psf_fit = synthpsf.(CartesianIndices(psf))
```

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
             loss=abs2,
             alg=NewtonTrustRegion(),
             maxfwhm=Inf,
             inv_var=nothing,
             kwargs...) where T
    _keys = keys(params)
    if isnothing(inv_var)
        @warn "Without inverse variance weights, cannot estimate PSF parameter covariance matrix." maxlog=5
    end
    _loss = build_loss_function(model, params, image, inv_var, inds; func_kwargs, loss, maxfwhm)
    X0 = vector_from_params(T, params)
    result = optimize(_loss, X0, alg, Optim.Options(; kwargs...); autodiff=ADTypes.AutoForwardDiff())
    Optim.converged(result) || @warn "optimizer did not converge" result
    X = Optim.minimizer(result)
    P_best = generate_params(_keys, X)
    best_model = model(T; P_best..., func_kwargs...)
    if isnothing(inv_var)
        return P_best, best_model
    end
    hessian = ForwardDiff.hessian(_loss, X)
    cov = inv(hessian)
    return P_best, best_model, cov
end

function build_loss_function(model::Model, params, image, inv_var, inds=axes(image); func_kwargs=(;), loss=abs2, maxfwhm=Inf)
    _keys = keys(params)
    cartinds = CartesianIndices(inds)
    minind = map(minimum, inds)
    maxind = map(maximum, inds)
    @inline function loss_function(X::AbstractVector{T}) where T
        P = generate_params(_keys, X)
        # position is within stamp
        minind[1] - 0.5 ≤ P.x ≤ maxind[1] + 0.5 || return T(Inf)
        minind[2] - 0.5 ≤ P.y ≤ maxind[2] + 0.5 || return T(Inf)
        # fwhm is non-negative and below max value
        if :fwhm in _keys
            all(0 .< P.fwhm .< maxfwhm) || return T(Inf)
        end
        # ratio is strictly (0, 1)
        if :ratio in _keys
            0 < P.ratio < 1 || return T(Inf)
        end
        # avoid circular degeneracy with theta
        if :theta in _keys
            -45 < P.theta < 45 || return T(Inf)
        end
        # sum of errors (by default, with L2 norm == chi-squared)
        chi2 = sum(cartinds) do idx
            resid = model(T, idx; P..., func_kwargs...) - image[idx]
            return isnothing(inv_var) ? loss(resid) : inv_var[idx] * loss(resid)
        end
        return chi2
    end
    return loss_function
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

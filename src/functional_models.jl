# Functional model interface

# factor for scaling radius in terms of the fwhm
const AIRY_PRE = π / 0.973

const BivariateLike = Union{<:Tuple, <:AbstractVector}

function rotate_point(dx, dy, theta)
    # generate rotation matrix
    # (theta is degrees CCW from x-axis)
    R = RotMatrix{2}(-deg2rad(theta))
    # rotate points
    return R * SA[dx, dy]
end

@doc raw"""
    gaussian([T=Float64], point; x, y, fwhm, amp=1, theta=0, bkg=0)
    gaussian([T=Float64], px, py; x, y, fwhm, amp=1, theta=0, bkg=0)

An unnormalized bivariate Gaussian distribution. The position can be specified
in `(x, y)` coordinates as a `Tuple`, `AbstractVector`, or as separate
arguments. If `theta` is given, the PSF will be rotated by `theta` degrees
counter-clockwise from the x-axis. If `bkg` is given it will be added as a
scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in
mind that `theta` has no effect for isotropic distributions and is degenerate
with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm`
tuple)

# Functional form
```math
f(x | x̂, \mathrm{FWHM}) = \exp[-4 \ln(2) ⋅ ||x - x̂|| / \mathrm{FWHM}^2]
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the
square-distance, and `FWHM` is the full width at half-maximum. If `FWHM` is a
scalar, the Gaussian distribution will be isotropic. If `FWHM` is a vector or
tuple, the weighting is applied along each axis (diagonal).
"""
function gaussian(T, px, py; x, y, fwhm, amp = 1, theta = 0, bkg = 0)
    flux = amp * (π * (fwhm isa BivariateLike ? prod(fwhm) : fwhm^2) / -GAUSS_PRE)
    model = if fwhm isa BivariateLike
        GaussianPSF(x, y, fwhm..., theta, flux, bkg)
    else
        !iszero(theta) && @warn "isotropic gaussian is not affected by non-zero rotation angle $theta"
        CircularGaussianPSF(x, y, fwhm, flux, bkg)
    end
    return convert(T, evaluate(model, px, py))
end

"""
    normal

An alias for [`gaussian`](@ref)
"""
const normal = gaussian


@doc raw"""
    airydisk([T=Float64], point; x, y, fwhm, ratio=0, amp=1, theta=0, bkg=0)
    airydisk([T=Float64], px, py; x, y, fwhm, ratio=0, amp=1, theta=0, bkg=0)

An unnormalized Airy disk. The position can be specified in `(x, y)` coordinates
as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta` is given,
the PSF will be rotated by `theta` degrees counter-clockwise from the x-axis. If
`bkg` is given it will be added as a scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in
mind that `theta` has no effect for isotropic distributions and is degenerate
with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm`
tuple)

If `ratio` is supplied, this will be the Airy pattern for a centrally-obscured
aperture (e.g., a Newtonian telescope). This has a slightly expanded functional
form, and in general the central Airy disk will be smaller and the first Airy
ring will be brighter.

# Functional form

The Airy disk is a distribution over the radius `r` (the square-Euclidean distance)

```math
f(x | x̂, \mathrm{FWHM}) = [ 2J₁(q) / q ]^2
```
where `J₁` is the first-order Bessel function of the first kind and
```math
q ≈ π r D / λ ≈ π r / (0.973 × \mathrm{FWHM})
```

If user a non-zero central obscuration via `ratio`, the functional form becomes

```math
f(x | x̂, \mathrm{FWHM}, ϵ) = [ 2J₁(q) / q - 2ϵJ₁(ϵq) / q ]^2 / (1 - ϵ^2)^2
```
where ``ϵ`` is the ratio (``0 ≤ ϵ < 1``).
"""
airydisk(T, px, py; x, y, fwhm, ratio = 0, amp = one(T), theta = 0, bkg = 0) =
    convert(T, _airydisk(px, py, x, y, fwhm, ratio, amp, theta, bkg))

# isotropic
function _airydisk(px, py, x, y, fwhm, ratio, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic airydisk is not affected by non-zero rotation angle $theta"
    # unnormalized airydisk
    r = sqrt(dx^2 + dy^2)
    # short-circuit
    iszero(r) && return amp + background
    q = AIRY_PRE * r / fwhm
    I1 = 2 * besselj1(q) / q
    iszero(ratio) && return amp * I1^2 + background
    I2 = 2 * ratio * besselj1(q * ratio) / q
    return amp * (I1 - I2)^2 + background
end

# bivariate
function _airydisk(px, py, x, y, fwhm::BivariateLike, ratio, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    if !iszero(theta)
        dx, dy = rotate_point(dx, dy, theta)
    end
    # unnormalized airy disk
    fwhmx, fwhmy = fwhm
    q = AIRY_PRE * sqrt((dx / fwhmx)^2 + (dy / fwhmy)^2)
    # short-circuit
    iszero(q) && return amp + background
    I1 = 2 * besselj1(q) / q
    iszero(ratio) && return amp * I1^2 + background
    I2 = 2 * ratio * besselj1(q * ratio) / q
    return amp * (I1 - I2)^2 + background
end


@doc raw"""
    moffat([T=Float64], point; x, y, fwhm, alpha=1, amp=1, theta=0, bkg=0)
    moffat([T=Float64], px, py; x, y, fwhm, alpha=1, amp=1, theta=0, bkg=0)

Two dimensional Moffat model. The position can be specified in `(x, y)`
coordinates as a `Tuple`, `AbstractVector`, or as separate arguments. If `theta`
is given, the PSF will be rotated by `theta` degrees counter-clockwise from the
x-axis. If `bkg` is given it will be added as a scalar to the PSF.

The `fwhm` can be a scalar (isotropic) or a vector/tuple (diagonal). Keep in
mind that `theta` has no effect for isotropic distributions and is degenerate
with the `fwhm` parameters (i.e., theta=90 is the same as reversing the `fwhm`
tuple)

# Functional form
```math
f(x | x̂, \mathrm{FWHM}, α) = A (1 + (||x - x̂|| / \mathrm{FWHM} / 2)^2)^{-α}
```
where `x̂` and `x` are position vectors (indices) `||⋅||` represents the
distance, and `FWHM` is the full width at half-maximum.
If `fwhm` is a vector or tuple, the weighting is applied along each axis.

Note that this function technically uses the half width at half-maximum, defined
as ``\mathrm{HWHM} = \mathrm{FWHM}/2``, but for compatibility with the other
models, `fwhm` is used as an input parameter instead.
"""
moffat(T, px, py; x, y, fwhm, alpha = 1, amp = one(T), theta = 0, bkg = 0) =
    convert(T, _moffat(px, py, x, y, fwhm, alpha, amp, theta, bkg))

# isotropic
function _moffat(px, py, x, y, fwhm, alpha, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    !iszero(theta) && @warn "isotropic moffat is not affected by non-zero rotation angle $theta"
    # unnormalized moffat
    gamma = _moffat_fwhm_to_gamma(fwhm, alpha)
    dist = (dx / gamma)^2 + (dy / gamma)^2
    return amp / (1 + dist)^alpha + background
end

# http://openafox.com/science/peak-function-derivations.html#moffat
_moffat_gamma_to_fwhm(gamma, alpha) = 2 * gamma * sqrt(2^(1 / alpha) - 1)
_moffat_fwhm_to_gamma(fwhm, alpha) = fwhm / (2 * sqrt(2^(1 / alpha) - 1))

# bivariate
function _moffat(px, py, x, y, fwhm::BivariateLike, alpha, amp, theta, background)
    # find offset from center
    dx = px - x
    dy = py - y
    # rotate
    if !iszero(theta)
        dx, dy = rotate_point(dx, dy, theta)
    end
    # unnormalized moffat
    gammax, gammay = _moffat_fwhm_to_gamma.(fwhm, alpha)
    # unnormalized moffat
    dist = (dx / gammax)^2 + (dy / gammay)^2
    return amp / (1 + dist)^alpha + background
end


# codegen for common functionality
# if you add a new model, make sure it gets added here
for model in (:gaussian, :airydisk, :moffat)
    @eval begin
        $model(px, py; kwargs...) = $model(Float64, px, py; kwargs...)
        $model(point::BivariateLike; kwargs...) = $model(point...; kwargs...)
        $model(T, point::BivariateLike; kwargs...) = $model(T, point...; kwargs...)
        $model(idx::CartesianIndex; kwargs...) = $model(idx.I; kwargs...)
        $model(T, idx::CartesianIndex; kwargs...) = $model(T, idx.I; kwargs...)
        $model(; kwargs...) = (point...) -> $model(_curried_point(point...); kwargs...)
        $model(::Type{T}; kwargs...) where {T} = (point...) -> $model(T, _curried_point(point...); kwargs...)
    end
end

_curried_point(P::BivariateLike) = P
_curried_point(point...) = Tuple(point)

############################################
# Fitting
############################################

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
- `params`: a named tuple of the parameters to fit and their default values,
  e.g., `(x=20, y=20, fwhm=3)`
- `image::AbstractMatrix{T}`: the data to fit the model to
- `inds`: the indices of the data to fit to, by default the entire array

# Keyword arguments
- `func_kwargs`: a named tuple of extra keyword arguments to pass to the model;
  this enables freezing a parameter, e.g., `(fwhm=3)`
- `loss`: the loss function to minimize, by default `abs2`, the chi-squared loss (L2 norm)
- `maxfwhm`: the maximum FWHM, can be a number or tuple
- `alg`: the optimization algorithm to use, by default `Optim.NewtonTrustRegion()`
- `inv_var`: inverse variance weights on the input `image`; if provided,
  must be the same size as `image`, and the covariance matrix of the best-fitting
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
function fit(
        model::Model,
        params,
        image::AbstractMatrix{T},
        inds = axes(image);
        func_kwargs = (;),
        loss = abs2,
        alg = Optim.NewtonTrustRegion(),
        maxfwhm = Inf,
        inv_var = nothing,
        kwargs...
    ) where {T}
    _keys = keys(params)
    if isnothing(inv_var)
        @warn "Without inverse variance weights, cannot estimate PSF parameter covariance matrix." maxlog = 5
    else
        size(inv_var) == size(image) || throw(ArgumentError("`inv_var` must be the same size as `image`"))
        any(<=(0), inv_var) && throw(ArgumentError("`inv_var` must be > 0 everywhere"))
    end
    _loss = build_loss_function(model, params, image, inv_var, inds; func_kwargs, loss, maxfwhm)
    X0 = vector_from_params(T, params)
    result = Optim.optimize(_loss, X0, alg, Optim.Options(; kwargs...); autodiff = ADTypes.AutoForwardDiff())
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

function build_loss_function(model::Model, params, image, inv_var, inds = axes(image); func_kwargs = (;), loss = abs2, maxfwhm = Inf)
    _keys = keys(params)
    cartinds = CartesianIndices(inds)
    minind = map(minimum, inds)
    maxind = map(maximum, inds)
    @inline function loss_function(X::AbstractVector{T}) where {T}
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
        first_half = @views _vals[begin:(_ind - 1)]
        fwhmx, fwhmy = _vals[_ind]
        second_half = @views _vals[(_ind + 1):end]
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
        first_half = @views zip(names[begin:(_ind - 1)], values[begin:(_ind - 1)])
        second_half = @views zip(names[(_ind + 1):end], values[(_ind + 2):end])
        return (; first_half..., fwhm, second_half...)
    end
    return (; zip(names, values)...)
end

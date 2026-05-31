
# ---------------------------------------------------------------------------
# Parameter vector ↔ model instance
# ---------------------------------------------------------------------------

"""
    free_params(model::AbstractPSFModel, fixed::NamedTuple=NamedTuple()) -> (free_names, x0)

Return the names of the free (non-fixed) parameters and their initial values as
a `Vector`. `fixed` is a `NamedTuple` whose keys name the fields to freeze.
"""
function free_params(model::AbstractPSFModel{T}, fixed::NamedTuple=NamedTuple()) where {T}
    all_props = getproperties(model)
    free_names = Tuple(k for k in keys(all_props) if !haskey(fixed, k))
    x0 = T[all_props[k] for k in free_names]
    return free_names, x0
end

"""
    model_from_vector(model, free_names, x, fixed) -> updated model

Reconstruct the model from an optimizer vector `x`, merging in the fixed
parameters.
"""
function model_from_vector(model, free_names::NTuple{N, Symbol}, x::AbstractVector, fixed::NamedTuple) where {N}
    updates = NamedTuple{free_names}(ntuple(i -> x[i], Val(N)))
    return setproperties(model, merge(updates, fixed))
end

# Type-stable internal variant: names are encoded in the Val type parameter so
# Julia can specialize NamedTuple{names} and setproperties concretely.
function model_from_vector(model, ::Val{names}, x::AbstractVector, fixed::NamedTuple) where {names}
    updates = NamedTuple{names}(ntuple(i -> x[i], Val(length(names))))
    return setproperties(model, merge(updates, fixed))
end

# ---------------------------------------------------------------------------
# Capability detection
# ---------------------------------------------------------------------------

function _has_hessian(model::AbstractPSFModel)
    return hasmethod(evaluate_fgh, Tuple{typeof(model), Real, Real})
end

function _has_deriv(model::AbstractPSFModel)
    return hasmethod(evaluate_fg, Tuple{typeof(model), Real, Real})
end

# ---------------------------------------------------------------------------
# NLSolversBase objective builders
# ---------------------------------------------------------------------------

function _make_fgh(model0::AbstractPSFModel{T}, free_names, free_idx, fixed, image, inds, inv_var, loss, maxfwhm) where {T}
    FT = float(T)
    n = length(free_idx)
    g_buf = Vector{FT}(undef, n)
    h_buf = Matrix{FT}(undef, n, n)
    free_names_val = Val(free_names)
    return NLSolversBase.only_fgh!(
        function (F, G, H, x)
            m = model_from_vector(model0, free_names_val, x, fixed)
            need_G = !isnothing(G)
            need_H = !isnothing(H)
            if need_G; fill!(G, 0); end
            if need_H; fill!(H, 0); end
            loss_val = zero(float(one(eltype(x))))
            for idx in CartesianIndices(inds)
                px, py = Tuple(idx)
                w = isnothing(inv_var) ? one(eltype(x)) : inv_var[idx]
                if need_H
                    f_val, g_full, h_full = evaluate_fgh(m, px, py)
                    d = image[idx]
                    ld  = LossFunctions.deriv(loss, f_val, d)
                    ldd = LossFunctions.deriv2(loss, f_val, d)
                    loss_val += w * loss(f_val, d)
                    @inbounds for (i, fi) in enumerate(free_idx)
                        g_buf[i] = g_full[fi]
                        for (j, fj) in enumerate(free_idx)
                            h_buf[i, j] = h_full[fi, fj]
                        end
                    end
                    wld  = w * ld
                    wldd = w * ldd
                    @inbounds for j in 1:n
                        gj = g_buf[j]
                        for i in 1:n
                            H[i, j] += wldd * g_buf[i] * gj + wld * h_buf[i, j]
                        end
                    end
                    if need_G
                        @inbounds for i in 1:n
                            G[i] += wld * g_buf[i]
                        end
                    end
                elseif need_G
                    f_val, g_full = evaluate_fg(m, px, py)
                    d = image[idx]
                    ld = LossFunctions.deriv(loss, f_val, d)
                    loss_val += w * loss(f_val, d)
                    wld = w * ld
                    @inbounds for (i, fi) in enumerate(free_idx)
                        G[i] += wld * g_full[fi]
                    end
                else
                    loss_val += w * loss(evaluate(m, px, py), image[idx])
                end
            end
            return !isnothing(F) ? loss_val : nothing
        end
    )
end

function _make_fg(model0::AbstractPSFModel{T}, free_names, free_idx, fixed, image, inds, inv_var, loss, maxfwhm) where {T}
    FT = float(T)
    n = length(free_idx)
    g_buf = Vector{FT}(undef, n)
    free_names_val = Val(free_names)
    return NLSolversBase.only_fg!(
        function (F, G, x)
            m = model_from_vector(model0, free_names_val, x, fixed)
            if !_bounds_ok(m, maxfwhm)
                if G !== nothing; fill!(G, 0); end
                return F !== nothing ? oftype(float(one(eltype(x))), Inf) : nothing
            end
            if G !== nothing; fill!(G, 0); end
            loss_val = zero(float(one(eltype(x))))
            for idx in CartesianIndices(inds)
                px, py = Tuple(idx)
                w = isnothing(inv_var) ? one(eltype(x)) : inv_var[idx]
                if G !== nothing
                    f_val, g_full = evaluate_fg(m, px, py)
                    d = image[idx]
                    ld = LossFunctions.deriv(loss, f_val, d)
                    loss_val += w * loss(f_val, d)
                    wld = w * ld
                    @inbounds for (i, fi) in enumerate(free_idx)
                        g_buf[i] = g_full[fi]
                    end
                    @inbounds for i in 1:n
                        G[i] += wld * g_buf[i]
                    end
                else
                    loss_val += w * loss(evaluate(m, px, py), image[idx])
                end
            end
            return F !== nothing ? loss_val : nothing
        end
    )
end

# ---------------------------------------------------------------------------
# Public fit method for struct-based models
# ---------------------------------------------------------------------------

"""
    PSFModels.fit(model::AbstractPSFModel, image, inds=axes(image);
                  fixed=(;), loss=LossFunctions.L2DistLoss(), maxfwhm=Inf,
                  alg=nothing, inv_var=nothing, kwargs...)

Fit the free parameters of `model` to `image[inds]` using `loss` as the
per-pixel loss function (default: L2 / chi-squared).

`fixed` is a `NamedTuple` of field-name → value pairs whose parameters are
frozen during the fit. All other fields of `model` are free.

Inverse variance weights can be passed via `inv_var`; if provided it must be
the same size as `image`. When `inv_var` is given the covariance matrix of the
free parameters is also returned.

Returns `(best_model, result)`, or `(best_model, result, cov)` when `inv_var`
is provided.

The optimizer is chosen automatically from the model's analytic-derivative
capabilities and the differentiability of `loss`:
- `evaluate_fgh` defined **and** `istwicedifferentiable(loss)` → `NewtonTrustRegion`
- `evaluate_fg` defined **and** `isdifferentiable(loss)` → `LBFGS`
- otherwise → `NelderMead`

Override by passing `alg` explicitly.
"""
function fit(model::AbstractPSFModel{T},
             image::AbstractMatrix,
             inds=axes(image);
             fixed::NamedTuple=(;),
             loss::LossFunctions.SupervisedLoss=LossFunctions.L2DistLoss(),
             maxfwhm=Inf,
             alg=nothing,
             inv_var=nothing,
             kwargs...) where {T}

    if !isnothing(inv_var)
        size(inv_var) == size(image) ||
            throw(ArgumentError("`inv_var` must be the same size as `image`"))
        any(<=(0), inv_var) &&
            throw(ArgumentError("`inv_var` must be > 0 everywhere"))
    end

    free_names, x0 = free_params(model, fixed)
    free_idx = findall(k -> !haskey(fixed, k), keys(getproperties(model)))

    use_hessian = _has_hessian(model) && LossFunctions.istwicedifferentiable(loss)
    use_deriv = _has_deriv(model) && LossFunctions.isdifferentiable(loss)

    alg = if !isnothing(alg)
        alg
    elseif use_hessian
        Optim.NewtonTrustRegion()
    elseif use_deriv
        Optim.LBFGS()
    else
        Optim.NelderMead()
    end

    if use_hessian || (!isnothing(alg) && _has_hessian(model))
        objective = _make_fgh(model, free_names, free_idx, fixed, image, inds, inv_var, loss, maxfwhm)
        result = Optim.optimize(objective, x0, alg, Optim.Options(; kwargs...))
    elseif use_deriv || (!isnothing(alg) && _has_deriv(model))
        objective = _make_fg(model, free_names, free_idx, fixed, image, inds, inv_var, loss, maxfwhm)
        result = Optim.optimize(objective, x0, alg, Optim.Options(; kwargs...))
    else
        # No analytic derivatives — scalar loss only (fall back to AD via Optim)
        function _scalar_loss(xv)
            m = model_from_vector(model, free_names, xv, fixed)
            _bounds_ok(m, maxfwhm) || return oftype(float(one(eltype(xv))), Inf)
            return sum(CartesianIndices(inds)) do idx
                px, py = Tuple(idx)
                w = isnothing(inv_var) ? one(eltype(xv)) : inv_var[idx]
                w * loss(evaluate(m, px, py), image[idx])
            end
        end
        result = Optim.optimize(_scalar_loss, x0, alg, Optim.Options(; kwargs...);
                                autodiff=ADTypes.AutoForwardDiff())
    end
    Optim.converged(result) || @warn "optimizer did not converge" result

    x_best = Optim.minimizer(result)
    best_model = model_from_vector(model, Val(free_names), x_best, fixed)

    if isnothing(inv_var)
        return best_model, result
    end

    # Covariance: H⁻¹ at the optimum, computed via ForwardDiff
    function _loss_for_hess(xv)
        m = model_from_vector(model, free_names, xv, fixed)
        return sum(CartesianIndices(inds)) do idx
            px, py = Tuple(idx)
            w = inv_var[idx]
            w * loss(evaluate(m, px, py), image[idx])
        end
    end
    hess = ForwardDiff.hessian(_loss_for_hess, x_best)
    cov = inv(hess)
    return best_model, result, cov
end

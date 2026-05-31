
# ---------------------------------------------------------------------------
# Parameter vector ↔ model instance
# ---------------------------------------------------------------------------

"""
    free_params(model::AbstractPSFModel, fixed::NamedTuple=NamedTuple()) -> (free_names, free_idx, x0)

Return the names of the free (non-fixed) parameters, their indices, and their initial values as
a `Vector`. `fixed` is a `NamedTuple` whose keys name the fields to freeze.
"""
function free_params(model::AbstractPSFModel{T}, fixed::NamedTuple=NamedTuple()) where {T}
    all_props = getproperties(model)
    free_names = Tuple(k for k in keys(all_props) if !haskey(fixed, k))
    free_idx   = Tuple(i for (i,k) in enumerate(keys(all_props)) if !haskey(fixed,k))
    x0 = T[all_props[k] for k in free_names]
    return free_names, free_idx, x0
end

"""
    model_from_vector(model, ::Val{names}, x, fixed::NamedTuple) -> updated model

Reconstruct the model from an optimizer vector `x`, merging in the fixed
parameters. `names` is a `Val`-wrapped tuple of the free parameter names
for type-stability, and `fixed` is a `NamedTuple` of the fixed parameters.
"""
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

function _make_fgh(model0::AbstractPSFModel{T}, free_names_val::Val, free_idx::Tuple, fixed,
                   image, inds, inv_var, loss) where {T}

    FT = float(T)

    return NLSolversBase.only_fgh!(
        function (F, G, H, x)
            m = model_from_vector(model0, free_names_val, x, fixed)
            if !isnothing(H)
                fill!(H, zero(eltype(H)))
                if !isnothing(G)
                    fill!(G, zero(eltype(G)))
                end
                loss_val = zero(FT)
                for idx in CartesianIndices(inds)
                    px = idx[1]
                    py = idx[2]
                    w = isnothing(inv_var) ? one(FT) : inv_var[idx]
                    f_val, g_full, h_full = evaluate_fgh(m, px, py)
                    d = image[idx]
                    ld  = LossFunctions.deriv(loss, f_val, d)
                    ldd = LossFunctions.deriv2(loss, f_val, d)
                    loss_val += w * loss(f_val, d)
                    wld  = w * ld
                    wldd = w * ldd
                    @inbounds for (j, fj) in pairs(free_idx)
                        gj = g_full[fj]
                        for (i, fi) in pairs(free_idx)
                            H[i, j] += wldd * g_full[fi] * gj +
                                       wld  * h_full[fi, fj]
                        end
                    end
                    if !isnothing(G)
                        @inbounds for (i, fi) in pairs(free_idx)
                            # G[i] += wld * g_full[fi]
                            G[i] = muladd(wld, g_full[fi], G[i])
                        end
                    end
                end
                return ifelse(!isnothing(F), loss_val, nothing)
            elseif !isnothing(G)
                fill!(G, zero(eltype(G)))
                loss_val = zero(FT)
                for idx in CartesianIndices(inds)
                    px = idx[1]
                    py = idx[2]
                    w = isnothing(inv_var) ? one(FT) : inv_var[idx]
                    f_val, g_full = evaluate_fg(m, px, py)
                    d = image[idx]
                    ld = LossFunctions.deriv(loss, f_val, d)
                    loss_val += w * loss(f_val, d)
                    wld = w * ld
                    @inbounds for (i, fi) in pairs(free_idx)
                        G[i] += wld * g_full[fi]
                    end
                end
                return ifelse(!isnothing(F), loss_val, nothing)
            else
                loss_val = zero(FT)
                for idx in CartesianIndices(inds)
                    px = idx[1]
                    py = idx[2]
                    w = isnothing(inv_var) ? one(FT) : inv_var[idx]
                    loss_val += w * loss(
                        evaluate(m, px, py),
                        image[idx]
                    )
                end
                return ifelse(!isnothing(F), loss_val, nothing)
            end
        end
    )
end

function _make_fg(model0::AbstractPSFModel{T}, free_names_val::Val, free_idx::Tuple, fixed,
                  image, inds, inv_var, loss) where {T}

    FT = float(T)

    return NLSolversBase.only_fg!(
        function (F, G, x)
            m = model_from_vector(model0, free_names_val, x, fixed)
            if !isnothing(G)
                fill!(G, zero(eltype(G)))
                loss_val = zero(FT)
                for idx in CartesianIndices(inds)
                    px = idx[1]
                    py = idx[2]
                    w = isnothing(inv_var) ? one(FT) : inv_var[idx]
                    f_val, g_full = evaluate_fg(m, px, py)
                    d = image[idx]
                    ld = LossFunctions.deriv(loss, f_val, d)
                    loss_val += w * loss(f_val, d)
                    wld = w * ld
                    @inbounds for (i, fi) in pairs(free_idx)
                        G[i] += wld * g_full[fi]
                    end
                end
                return !isnothing(F) ? loss_val : nothing
            else
                loss_val = zero(FT)
                for idx in CartesianIndices(inds)
                    px = idx[1]
                    py = idx[2]
                    w = isnothing(inv_var) ? one(FT) : inv_var[idx]
                    loss_val += w * loss(
                        evaluate(m, px, py),
                        image[idx],
                    )
                end
                return ifelse(!isnothing(F), loss_val, nothing)
            end
        end,
    )
end

# ---------------------------------------------------------------------------
# Public fit method for struct-based models
# ---------------------------------------------------------------------------

"""
    PSFModels.fit(model::AbstractPSFModel, image, inds=axes(image);
                  fixed=(;), loss=LossFunctions.L2DistLoss(),
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
             alg=nothing,
             inv_var=nothing,
             kwargs...) where {T}

    if !isnothing(inv_var)
        size(inv_var) == size(image) ||
            throw(ArgumentError("`inv_var` must be the same size as `image`"))
        any(<=(0), inv_var) &&
            throw(ArgumentError("`inv_var` must be > 0 everywhere"))
    end

    free_names, free_idx, x0 = free_params(model, fixed)
    free_names_val = Val(free_names) # for type stability in model_from_vector

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

    objective = if alg isa Union{Optim.FirstOrderOptimizer, Optim.SecondOrderOptimizer} && use_hessian
                    _make_fgh(model, free_names_val, free_idx, fixed, image, inds, inv_var, loss)
                elseif alg isa Optim.FirstOrderOptimizer && use_deriv
                    _make_fg(model, free_names_val, free_idx, fixed, image, inds, inv_var, loss)
                else
                    # No analytic derivatives — scalar loss only (fall back to AD via Optim)
                    function _scalar_loss(xv)
                        m = model_from_vector(model, free_names_val, xv, fixed)
                        return sum(CartesianIndices(inds)) do idx
                            px, py = idx[1], idx[2]
                            w = isnothing(inv_var) ? one(eltype(xv)) : inv_var[idx]
                            w * loss(evaluate(m, px, py), image[idx])
                        end
                    end
                end
    result = Optim.optimize(objective, x0, alg, Optim.Options(; kwargs...); autodiff=ADTypes.AutoForwardDiff())
    Optim.converged(result) || @warn "optimizer did not converge" result

    x_best = Optim.minimizer(result)
    best_model = model_from_vector(model, free_names_val, x_best, fixed)

    if isnothing(inv_var)
        return best_model, result
    end

    # Covariance: H⁻¹ at the optimum, computed via ForwardDiff
    # TODO: if the model has an analytic Hessian, use it instead of ForwardDiff to save time and improve precision
    function _loss_for_hess(xv)
        m = model_from_vector(model, free_names_val, xv, fixed)
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


"""
    LMResult{T}

Result from [`PSFModels.fit_lm`](@ref).

Fields:
- `minimizer`:   final free-parameter vector
- `minimum`:     final cost `∑ wᵢ (fᵢ − dᵢ)²` at the solution
- `cost_init`:   cost at the initial parameter vector
- `converged`:   `true` when any termination criterion was satisfied
- `x_converged`: `true` when the step norm fell below `x_tol`
- `f_converged`: `true` when the cost decrease fell below `f_tol`
- `g_converged`: `true` when the gradient norm fell below `g_tol`
- `iterations`:  total number of Levenberg-Marquardt iterations performed
- `λ_final`:     damping parameter value at termination
- `σ_final`:     final scale estimate (NaN if not applicable)
- `cov`:         covariance matrix of the free parameters
"""
struct LMResult{T,V<:AbstractVector{T},M<:AbstractMatrix{T}}
    minimizer::V
    minimum::T
    cost_init::T
    converged::Bool
    x_converged::Bool
    f_converged::Bool
    g_converged::Bool
    iterations::Int
    λ_final::T
    σ_final::T
    cov::M
end

function Base.show(io::IO, r::LMResult)
    println(io, "LMResult:")
    println(io, "  converged:  ", r.converged,
                "  (x: ", r.x_converged,
                ", f: ", r.f_converged,
                ", g: ", r.g_converged, ")")
    println(io, "  iterations: ", r.iterations)
    println(io, "  cost:       ", r.minimum, "  (init: ", r.cost_init, ")")
    if !isnan(r.σ_final)
        println(io, "  σ_final:    ", r.σ_final)
    end
    print(io,   "  λ_final:    ", r.λ_final)
end

# ---------------------------------------------------------------------------
# LM damping strategies
# ---------------------------------------------------------------------------

abstract type AbstractLMDamping end
"""
    MarquardtDamping(; min_diagonal=1e-6)::MarquardtDamping

Damping strategy for Levenberg-Marquardt where the diagonal entries are
scaled by `max(A[i, i], min_diagonal)` to prevent small curvature directions from being under-damped.
"""
Base.@kwdef struct MarquardtDamping{T} <: AbstractLMDamping
    min_diagonal::T = 1e-6
end
function damp!(A, damping::MarquardtDamping, λ)
    min_diagonal = damping.min_diagonal
    @inbounds for i in 1:size(A, 1)
        A[i, i] += λ * max(A[i, i], min_diagonal)
    end
end

"""
    LevenbergDamping()::LevenbergDamping

Damping strategy for Levenberg-Marquardt where a uniform `λ I` shift is applied with no scaling.
"""
struct LevenbergDamping{T} <: AbstractLMDamping end
function damp!(A, damping::LevenbergDamping, λ)
    @inbounds for i in 1:size(A, 1)
        A[i, i] += λ
    end
end

"""
    NoDamping()::NoDamping

Damping strategy for Levenberg-Marquardt where no damping is applied; equivalent to Gauss-Newton. Not recommended for general use.
"""
struct NoDamping <: AbstractLMDamping end
function damp!(A, damping::NoDamping, λ) end # No-op

# ---------------------------------------------------------------------------
# Scale estimators for IRLS
# ---------------------------------------------------------------------------

"""Abstract type for robust scale estimators used in iteratively reweighted least squares."""
abstract type AbstractScaleEstimator end

"""
    MADScale()

Scale estimator using the median absolute deviation:
``σ̂ = 1.4826 · median(|rᵢ − median(rᵢ)|)``

The factor 1.4826 makes the estimate consistent for normally distributed data.

# Examples

```jldoctest
julia> using PSFModels: MADScale, estimate_scale

julia> r = randn(10_000);

julia> isapprox(estimate_scale(MADScale(), r), 1.0; atol=0.05)
true

julia> r2 = [randn(9_000); 100.0 * ones(1_000)]  # 10% outliers

julia> isapprox(estimate_scale(MADScale(), r2), 1.15; atol=0.05)
true
```
"""
struct MADScale <: AbstractScaleEstimator end

function estimate_scale(::MADScale, r::AbstractVector{<:Real})
    med = median(r)
    return 1.4826 * median(abs.(r .- med))
end

"""
    FixedScale(σ::Real)

Scale estimator that returns a fixed, user-provided scale value.

# Examples

```jldoctest
julia> using PSFModels: FixedScale, estimate_scale

julia> estimate_scale(FixedScale(2.0), randn(100))
2.0
```
"""
struct FixedScale{T<:Real} <: AbstractScaleEstimator
    σ::T
end
estimate_scale(est::FixedScale, ::AbstractVector) = est.σ

"""
    MScale(; δ=0.5, tol=1e-6, max_iter=30)

Iterative M-scale estimator solving ``(1/n) ∑ χ(rᵢ/σ) = δ`` where
``χ(r) = r²`` (Huber's proposal 2). The breakdown point is controlled
by ``δ`` (default 0.5 for 50% breakdown). Iteration stops when the
relative change in ``σ`` falls below ``tol`` or ``max_iter`` is reached.

For Gaussian data with ``δ = 0.5``, the estimator returns ``σ ≈ √2 = 1.41``
since ``E[χ(r/σ)] ≈ 1/σ²`` and the equation requires this to equal ``δ``.
This makes ``MScale`` consistent with the IRLS weight-function thresholds,
which are calibrated to the scale *after* M-estimation.

# Examples

```jldoctest
julia> using PSFModels: MScale, estimate_scale

julia> r = randn(1_000);

julia> isapprox(estimate_scale(MScale(), r), sqrt(2); atol=0.1)
true
```
"""
Base.@kwdef struct MScale{T} <: AbstractScaleEstimator
    δ::T = 0.5
    tol::T = 1e-6
    max_iter::Int = 30
end

function estimate_scale(est::MScale, r::AbstractVector{<:Real})
    T = float(eltype(r))
    δ = T(est.δ)
    tol = T(est.tol)
    # Initial estimate from MAD
    σ = estimate_scale(MADScale(), r)
    n = length(r)
    for _ in 1:est.max_iter
        σ2 = σ^2
        # χ(r) = r² for |r| ≤ 3σ, else (3σ)² (capped to bound influence)
        cap = 9 * σ2
        χ_sum = zero(T)
        @inbounds for ri in r
            ri2 = ri^2
            χ_sum += ifelse(ri2 < cap, ri2, cap)
        end
        σ_new = sqrt(χ_sum / (n * δ))
        if abs(σ_new - σ) ≤ tol * max(σ, tol)
            σ = σ_new
            break
        end
        σ = σ_new
    end
    return σ
end

# ---------------------------------------------------------------------------
# Custom loss types for IRLS
# ---------------------------------------------------------------------------

"""
    TukeyLoss(; c=4.685)

Tukey's bisquare (biweight) loss function for robust estimation.

The loss is bounded, meaning large residuals are completely rejected
(weight → 0 for |r| ≥ c).  This makes it ideal for rejecting cosmic rays
and other severe outliers in astronomical images.

The default tuning constant ``c = 4.685`` gives 95% asymptotic efficiency
under Gaussian errors.  Lower values provide more aggressive outlier
rejection at the cost of Gaussian efficiency.

# Examples

```jldoctest
julia> using PSFModels: TukeyLoss

julia> using LossFunctions

julia> loss = TukeyLoss(; c=4.685);

julia> loss(0.0, 0.0)
0.0

julia> loss(3.0, 1.0) > 0
true

julia> isapprox(LossFunctions.deriv(loss, 0.1, 0.0), 0.1; atol=1e-3) # for |r| << c, ψ(r) ≈ r
true

julia> isapprox(LossFunctions.deriv(loss, 0.3, 0.1), 0.2; atol=1e-3) # for |r| << c, ψ(r) ≈ r
true

julia> ψ_expected = 2.0 * (1 - (2.0 / 4.685)^2)^2;

julia> isapprox(LossFunctions.deriv(loss, 3.0, 1.0), ψ_expected; atol=1e-6)
true

julia> isapprox(LossFunctions.deriv(loss, 10.0, 0.0), 0.0; atol=1e-6) # beyond threshold
true

julia> isapprox(LossFunctions.deriv2(loss, 0.0, 0.0), 1.0; atol=1e-6) # Ψ'(0) = 1
true

julia> LossFunctions.deriv2(loss, 10.0, 0.0) == 0.0 # Ψ'(r) = 0 for |r| ≥ c
true
```
"""
Base.@kwdef struct TukeyLoss{T} <: LossFunctions.SupervisedLoss
    c::T = 4.685
end

# ρ(r) = (c²/6) · [1 − (1 − (r/c)²)³]  for |r| ≤ c,  else c²/6
function (loss::TukeyLoss)(y, t)
    T = promote_type(typeof(y), typeof(t))
    T = float(T)
    r = y - t
    c = T(loss.c)
    if abs(r) ≤ c
        r_over_c2 = (r / c)^2
        return (c^2 / 6) * (1 - (1 - r_over_c2)^3)
    else
        return c^2 / 6
    end
end

# ψ(r) = r · (1 − (r/c)²)²  for |r| ≤ c,  else 0
function LossFunctions.deriv(loss::TukeyLoss, y, t)
    T = promote_type(typeof(y), typeof(t))
    T = float(T)
    c = T(loss.c)
    r = y - t
    if abs(r) ≤ c
        r_over_c2 = (r / c)^2
        return r * (1 - r_over_c2)^2
    else
        return zero(T)
    end
end

# ψ'(r) = (1 − (r/c)²) · (1 − 5(r/c)²)  for |r| ≤ c,  else 0
function LossFunctions.deriv2(loss::TukeyLoss, y, t)
    T = promote_type(typeof(y), typeof(t))
    T = float(T)
    c = T(loss.c)
    r = y - t
    if abs(r) ≤ c
        r_over_c2 = (r / c)^2
        return (1 - r_over_c2) * (1 - 5 * r_over_c2)
    else
        return zero(T)
    end
end

# ---------------------------------------------------------------------------
# IRLS weight helper
# ---------------------------------------------------------------------------

"""
    weight(loss::LossFunctions.SupervisedLoss, r::Real) → Real

Compute the IRLS weight ``w(r) = ψ(r) / (r · ψ'(0))`` where ``ψ`` is the
influence function (first derivative) of `loss`.  The normalization by
``ψ'(0)`` ensures that ``w(0) = 1`` for any loss function.

# Examples

```jldoctest
julia> using LossFunctions: L2DistLoss, HuberLoss

julia> using PSFModels: TukeyLoss, weight

julia> weight(L2DistLoss(), 1.0) == 1.0
true

julia> weight(HuberLoss(1.0), 0.5) == 1.0  # within threshold
true

julia> weight(TukeyLoss(; c=4.685), 10.0) == 0.0  # beyond threshold
true
```
"""
function weight(loss::LossFunctions.SupervisedLoss, r::T) where {T}
    ψ0 = LossFunctions.deriv2(loss, zero(T), zero(T))
    if abs(r) < eps(T) * 10
        return one(T)
    end
    return LossFunctions.deriv(loss, r, zero(T)) / (r * ψ0)
end

# ---------------------------------------------------------------------------
# fit_lm
# ---------------------------------------------------------------------------

# Build the Gauss-Newton approximation to the Hessian (A = Jᵀ W J) and the
# gradient (b = Jᵀ W r) together with the cost C = ∑ wᵢ rᵢ².
#
# The Jacobian row for pixel i is `g_full[free_idx]`, the projection of the
# full analytic gradient onto the free parameters.  g_full is a StaticArrays
# SVector (immutable), so we index into it rather than mutating it.
function _lm_accum_fg!(A::AbstractMatrix{FT}, b::AbstractVector{FT}, m::AbstractPSFModel, free_idx, image, inds, inv_var, weights=nothing) where {FT}
    @assert size(A, 1) == size(A, 2) == length(b) == length(free_idx)
    fill!(A, zero(FT))
    fill!(b, zero(FT))
    cost = zero(FT)
    n = length(free_idx)
    use_weights = weights !== nothing
    for idx in CartesianIndices(inds)
        px, py = idx[1], idx[2]
        w = if use_weights
            FT(weights[idx])
        elseif isnothing(inv_var)
            one(FT)
        else
            FT(inv_var[idx])
        end
        f_val, g_full = evaluate_fg(m, px, py)
        r  = FT(f_val) - FT(image[idx])
        wr = w * r
        cost = muladd(w * r, r, cost)
        @inbounds for j in 1:n
            gj = FT(g_full[free_idx[j]])
            b[j] = muladd(wr, gj, b[j])
            for i in 1:n
                A[i, j] = muladd(w * FT(g_full[free_idx[i]]), gj, A[i, j])
            end
        end
    end
    return cost
end

"""
    PSFModels.fit_lm(model::AbstractPSFModel, image, inds=axes(image);
                     fixed=(;), inv_var=nothing,
                     damping::AbstractLMDamping=MarquardtDamping(),
                     max_iter=200, x_tol=1e-8, f_tol=1e-8, g_tol=1e-8,
                     show_trace=false, reweight=nothing,
                     scale_estimator=nothing)

Fit the free parameters of `model` to `image[inds]` under weighted L2 loss
using the Levenberg-Marquardt algorithm.  The model must implement
`evaluate_fg`.

`fixed` is a `NamedTuple` of field-name → value pairs whose parameters are
frozen during the fit.  All other fields of `model` are free.

Inverse variance weights can be passed via `inv_var`; it must be the same size
as `image`.

Returns an `LMResult` containing the fitted parameters, cost, convergence
flags, and covariance matrix.

# Algorithm

Minimises ``C(\\mathbf{x}) = \\sum_i w_i [f_i(\\mathbf{x}) - d_i]^2`` using
the Gauss-Newton approximation to the Hessian augmented by a diagonal damping
term:

```math
(J^\\top W J + \\lambda D)\\,\\delta = -J^\\top W r
```

where ``J`` is the Jacobian assembled from `evaluate_fg` gradients (restricted
to the free parameters), ``W = \\mathrm{diag}(w_i)``, and
``r_i = f_i(\\mathbf{x}) - d_i``.  If the candidate step reduces the cost it
is accepted and ``\\lambda`` is decreased; otherwise it is rejected and
``\\lambda`` is increased.

# Iteratively Reweighted Least Squares

Pass `reweight` as a `LossFunctions.SupervisedLoss` (e.g. `HuberLoss()`,
`TukeyLoss()`) to enable iteratively reweighted least squares (IRLS). After
each accepted LM step, the residuals are used to recompute pixel weights via
``w_i^{\\text{final}} = w_i^{\\text{base}} \\cdot w(r_i/\\sigma)`` where
``w(r) = \\psi(r)/r`` is derived from the loss function's influence function
and ``\\sigma`` is a robust scale estimate. This provides automatic outlier
rejection — useful for cosmic rays, bad pixels, and satellite trails in
astronomical images.

The scale ``\\sigma`` is estimated by `scale_estimator`. If not provided,
it defaults to `FixedScale` (from `inv_var`) when `inv_var` is supplied,
otherwise `MADScale()`. Available scale estimators:
- `MADScale()` — median absolute deviation (robust, default without `inv_var`)
- `FixedScale(σ)` — fixed user-provided scale
- `MScale(; δ=0.5)` — iterative M-scale (most robust, highest cost)

Recommended loss functions for IRLS (all from `LossFunctions` or `PSFModels`):
- `HuberLoss(c)` — soft downweighting of outliers, suitable for mildly
  contaminated data or crowded-field photometry
- `TukeyLoss(; c=4.685)` — complete rejection of extreme outliers, ideal
  for cosmic rays and bad pixels in space-based imaging

# Keyword arguments

- `fixed`: `NamedTuple` of frozen parameter name → value pairs
- `inv_var`: inverse-variance weights, same shape as `image`
- `reweight`: a `LossFunctions.SupervisedLoss` for IRLS reweighting,
  or `nothing` (default) for standard L2 fitting
- `scale_estimator::AbstractScaleEstimator`: scale estimator for IRLS;
  defaults to `FixedScale` (from `inv_var`) or `MADScale()` (see above)
- `λ_init`: initial damping parameter (default `1e-4`)
- `λ_up`: factor by which `λ` is multiplied on rejection (default `10`)
- `λ_down`: factor by which `λ` is divided on acceptance (default `10`)
- `λ_min`, `λ_max`: lower/upper bounds on the damping parameter
- `damping::AbstractLMDamping`: specifies the damping strategy
  (`MarquardtDamping`, `LevenbergDamping`, or `NoDamping`)
- `max_iter`: maximum number of LM iterations (default `200`)
- `x_tol`: step-norm convergence criterion:
  ``\\|\\delta\\| \\le x\\_tol\\cdot(\\|x\\|+x\\_tol)`` (default `1e-8`)
- `f_tol`: cost-decrease convergence criterion:
  ``\\Delta C \\le f\\_tol\\cdot(|C|+f\\_tol)`` (default `1e-8`)
- `g_tol`: gradient-norm convergence criterion ``\\|J^\\top W r\\|``
  (default `1e-8`)
- `show_trace`: print per-iteration statistics to stdout (default `false`)

# Damping strategies
- `Marquardt`: diagonal entries are scaled by `max(A[i, i], min_diagonal)` to
  prevent small curvature directions from being under-damped
- `Levenberg`: uniform damping with no scaling
- `NoDamping`: no damping; equivalent to Gauss-Newton (not recommended)
"""
function fit_lm(model::AbstractPSFModel{T},
                image::AbstractMatrix,
                inds=axes(image);
                fixed::NamedTuple=(;),
                inv_var=nothing,
                λ_init::Real=1e-4,
                λ_up::Real=10.0,
                λ_down::Real=10.0,
                λ_min::Real=1e-12,
                λ_max::Real=1e12,
                damping::AbstractLMDamping=MarquardtDamping(),
                max_iter::Integer=200,
                x_tol::Real=1e-8,
                f_tol::Real=1e-8,
                g_tol::Real=1e-8,
                show_trace::Bool=false,
                reweight::Union{Nothing, LossFunctions.SupervisedLoss}=nothing,
                scale_estimator::Union{Nothing, AbstractScaleEstimator}=nothing) where {T}

    # Validate inputs
    if !isnothing(inv_var)
        size(inv_var) == size(image) ||
            throw(ArgumentError("`inv_var` must be the same size as `image`"))
        any(<=(0), inv_var) &&
            throw(ArgumentError("`inv_var` must be > 0 everywhere"))
    end
    _has_deriv(model) ||
        throw(ArgumentError("model does not implement `evaluate_fg`; " *
                            "Levenberg-Marquardt requires gradient"))

    # Parameter bookkeeping
    free_names, free_idx, x0 = free_params(model, fixed)
    n = length(x0)
    n > 0 || throw(ArgumentError("all model parameters are fixed; nothing to fit"))
    free_names_val = Val(free_names)

    # Pre-allocate for in-place accumulation
    FT = float(T)
    λ = FT(λ_init)
    x, x_cand = Vector{FT}(x0), Vector{FT}(undef, n)
    A, b = zeros(FT, n, n), zeros(FT, n)
    A_cand, b_cand = zeros(FT, n, n), zeros(FT, n)
    A_damp = zeros(FT, n, n)
    F_buffer = zeros(FT, n, n)
    δ = zeros(FT, n)

    # IRLS setup
    do_irls = reweight !== nothing
    # Determine scale estimator: explicit > FixedScale from inv_var > MAD
    scale_est = if !isnothing(scale_estimator)
        scale_estimator
    elseif !isnothing(inv_var)
        # FixedScale from inv_var — σ ≈ 1/√(inv_var) averaged over the cutout
        FixedScale(sqrt(FT(1) / mean(x -> FT(x), view(inv_var, inds...))))
    else
        MADScale()
    end
    σ_final = FT(NaN)
    # Pre-allocate weight arrays for IRLS
    weights = do_irls ? similar(image, FT) : nothing
    base_weights = do_irls ? similar(image, FT) : nothing
    residuals = do_irls ? Vector{FT}(undef, prod(length, inds)) : nothing
    if do_irls
        # Initialize base weights from inv_var (or 1)
        for idx in CartesianIndices(inds)
            base_weights[idx] = isnothing(inv_var) ? one(FT) : FT(inv_var[idx])
        end
        weights .= base_weights
    end

    # Initial cost and Jacobian
    m0 = model_from_vector(model, free_names_val, x, fixed)
    cost = _lm_accum_fg!(A, b, m0, free_idx, image, inds, inv_var, weights)
    cost_init  = cost

    if show_trace
        gnorm0 = sqrt(sum(abs2, b))
        println("Initialization | cost = $cost_init | λ = $λ | ||g|| = $gnorm0")
    end

    # LM iteration
    x_converged = false
    f_converged = false
    g_converged = false
    converged   = false
    iter = 0

    while iter < max_iter
        iter += 1

        # Test for gradient-norm convergence
        gnorm = sqrt(sum(abs2, b))
        if gnorm ≤ FT(g_tol)
            g_converged = true
            converged   = true
            show_trace && println("Iter $(lpad(iter,4)) | converged on gradient norm (||g|| = $gnorm)")
            break
        end

        # Build the damped normal matrix A + λD
        A_damp .= A
        damp!(A_damp, damping, λ)

        # Solve (A + λD) δ = −b; same as δ = A_damp \ (-b)
        F_buffer .= A_damp
        F = cholesky!(Symmetric(F_buffer))
        ldiv!(δ, F, -b)
        δnorm = sqrt(sum(abs2, δ))

        # Evaluate the candidate step
        x_cand .= x .+ δ
        m_cand = model_from_vector(model, free_names_val, x_cand, fixed)
        cost_cand = _lm_accum_fg!(A_cand, b_cand, m_cand, free_idx, image, inds, inv_var, weights)
        accepted = cost_cand < cost

        if show_trace
            status = accepted ? "accepted" : "rejected"
            println("Iter $(lpad(iter,4)) | cost = $cost → $cost_cand | " *
                    "λ = $λ | ||g|| = $gnorm | ||δ|| = $δnorm | $status")
        end

        if accepted
            Δcost = cost - cost_cand
            x    .= x_cand
            cost  = cost_cand
            A    .= A_cand
            b    .= b_cand
            λ    = max(λ / FT(λ_down), FT(λ_min))

            # Check standard LM convergence before reweighting
            x_converged = δnorm ≤ FT(x_tol) * (sqrt(sum(abs2, x)) + FT(x_tol))
            f_converged = Δcost ≤ FT(f_tol) * (abs(cost) + FT(f_tol))

            if x_converged || f_converged
                converged = true
                break
            end

            # IRLS: recompute weights from current residuals
            if do_irls
                # Compute residuals for the current model
                ri = 1
                for idx in CartesianIndices(inds)
                    px, py = idx[1], idx[2]
                    f_val = evaluate(m_cand, px, py)
                    residuals[ri] = FT(f_val) - FT(image[idx])
                    ri += 1
                end
                σ_final = estimate_scale(scale_est, residuals)

                # Skip reweighting if scale is zero (noiseless or degenerate)
                if σ_final > eps(FT)
                    ri = 1
                    for idx in CartesianIndices(inds)
                        r_scaled = residuals[ri] / σ_final
                        weights[idx] = base_weights[idx] * weight(reweight, FT(r_scaled))
                        ri += 1
                    end
                    # Recompute cost and gradient with new weights
                    cost = _lm_accum_fg!(A, b, m_cand, free_idx, image, inds, inv_var, weights)
                    # Reset λ — the Gauss-Newton approximation changed with the new weights
                    λ = FT(λ_init)
                end
            end
        else
            λ = min(λ * FT(λ_up), FT(λ_max))
        end
    end

    best_model = model_from_vector(model, free_names_val, x, fixed)

    # Compute covariance matrix
    function _compute_cov(JTJ, cost_val, n_obs, n_params)
        dof = n_obs - n_params
        dof > 0 || throw(ArgumentError("degrees of freedom must be positive to compute covariance"))
        σ² = cost_val / dof
        cov = try
            cholesky!(Symmetric(JTJ)) \ I
        catch
            # Fall back to pinv if Cholesky fails (may happen with heavy IRLS downweighting)
            pinv(Symmetric(JTJ))
        end
        cov *= σ²
        return cov
    end
    n_obs = prod(length, inds)
    cov = _compute_cov(A, cost, n_obs, n)

    result = LMResult(x, cost, cost_init,
                      converged, x_converged, f_converged, g_converged,
                      iter, λ, σ_final, cov)
    return result
end

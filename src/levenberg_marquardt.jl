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
- `chisq`:       final reduced chi-squared (cost per degree of freedom)
"""
struct LMResult{T, V <: AbstractVector{T}, M <: AbstractMatrix{T}}
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
    chisq::T
end

function Base.show(io::IO, r::LMResult)
    println(io, "LMResult:")
    println(
        io, "  converged:  ", r.converged,
        "  (x: ", r.x_converged,
        ", f: ", r.f_converged,
        ", g: ", r.g_converged, ")"
    )
    println(io, "  iterations: ", r.iterations)
    println(io, "  cost:       ", r.minimum, "  (init: ", r.cost_init, ")")
    println(io, "  reduced χ²:  ", r.chisq)
    if !isnan(r.σ_final)
        println(io, "  σ_final:    ", r.σ_final)
    end
    return print(io, "  λ_final:    ", r.λ_final)
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
    min_diagonal::T = 1.0e-6
end
function damp!(A, damping::MarquardtDamping, λ)
    min_diagonal = damping.min_diagonal
    return @inbounds for i in axes(A, 1)
        A[i, i] += λ * max(A[i, i], min_diagonal)
    end
end

"""
    LevenbergDamping()::LevenbergDamping

Damping strategy for Levenberg-Marquardt where a uniform `λ I` shift is applied with no scaling.
"""
struct LevenbergDamping <: AbstractLMDamping end
function damp!(A, damping::LevenbergDamping, λ)
    return @inbounds for i in axes(A, 1)
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

function estimate_scale(::MADScale, r::AbstractArray)
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
struct FixedScale{T <: Real} <: AbstractScaleEstimator
    σ::T
end
estimate_scale(est::FixedScale, ::AbstractArray) = est.σ

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
    tol::T = 1.0e-6
    max_iter::Int = 30
end

function estimate_scale(est::MScale, r::AbstractArray{T}) where {T}
    FT = float(T)
    δ = FT(est.δ)
    tol = FT(est.tol)
    # Initial estimate from MAD
    σ = max(estimate_scale(MADScale(), r), eps(FT))
    n = length(r)
    for _ in 1:est.max_iter
        σ2 = σ^2
        # χ(r) = r² for |r| ≤ 3σ, else (3σ)² (capped to bound influence)
        cap = 9 * σ2
        χ_sum = zero(FT)
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
# Covariance estimators
# ---------------------------------------------------------------------------

abstract type AbstractCovarianceEstimator end

"""`KnownWeightsCovarianceEstimator()` assumes that the weights provided (e.g. via `inv_var`) are correct and returns the covariance as the inverse of the Gauss-Newton Hessian approximation.
"""
struct KnownWeightsCovarianceEstimator <: AbstractCovarianceEstimator end
function covariance!(::KnownWeightsCovarianceEstimator, JTJ, cost_val, dof)
    # For known weights (e.g. from inv_var), the covariance
    # is simply the inverse of the Gauss-Newton Hessian approximation
    cov = try
        F = cholesky!(Symmetric(JTJ))
        F \ I # = inv(JTJ), more stable
    catch
        pinv(JTJ) # fallback to pseudo-inverse if JTJ is not positive definite
    end
    return cov
end
"""`ReweightedCovarianceEstimator()` inflates the covariance by the reduced cost per degree of freedom to account for the fact that the IRLS weights are estimated from the data and may not be correct."""
struct ReweightedCovarianceEstimator <: AbstractCovarianceEstimator end
function covariance!(::ReweightedCovarianceEstimator, JTJ, cost_val, dof)
    # For reweighted estimates, the covariance is inflated by the reduced cost per degree of freedom
    # F = cholesky!(Symmetric(JTJ))
    # cov = F \ I # = inv(JTJ), more stable
    cov = try
        F = cholesky!(Symmetric(JTJ))
        F \ I # = inv(JTJ), more stable
    catch
        pinv(JTJ) # fallback to pseudo-inverse if JTJ is not positive definite
    end
    return (cost_val / dof) * cov
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
function _lm_accum_fg!(A::AbstractMatrix{FT}, b::AbstractVector{FT}, m::AbstractPSFModel, free_idx, image, inds, residuals, inv_var, weights = nothing) where {FT}
    # Get the weight based on arguments provided
    function _w(use_weights, use_inv_var, weights, inv_var, idx, ::Type{FT}) where {FT}
        return FT(use_weights ? weights[idx] : (use_inv_var ? inv_var[idx] : one(FT)))
    end
    @assert size(A, 1) == size(A, 2) == length(b) == length(free_idx)
    @assert size(image) == size(residuals)
    fill!(A, zero(FT))
    fill!(b, zero(FT))
    cost = zero(FT)
    n = length(free_idx)
    # Weights come from IRLS, supercede inv_var if provided
    use_weights = !isnothing(weights)
    use_inv_var = !isnothing(inv_var)
    for idx in CartesianIndices(inds)
        w = _w(use_weights, use_inv_var, weights, inv_var, idx, FT)
        f_val, g_full = evaluate_fg(m, idx)
        r = FT(f_val) - FT(image[idx])
        residuals[idx] = r
        wr = w * r
        cost = muladd(wr, r, cost)
        @inbounds for j in 1:n
            gj = FT(g_full[free_idx[j]])
            b[j] = muladd(wr, gj, b[j])
            # A is symmetric so we could only compute upper triangle,
            # but this is not a bottleneck
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
                     show_trace=false, 
                     reweight=nothing,
                     scale_estimator=nothing,
                     covariance_estimator=nothing)

Fit the free parameters of `model` to `image[inds]` under weighted L2 loss
using the Levenberg-Marquardt algorithm.  The model must implement
`evaluate_fg`.

`fixed` is a `NamedTuple` of field-name → value pairs whose parameters are
frozen during the fit.  All other fields of `model` are free.

Inverse variance weights can be passed via `inv_var`; it must be the same size
as `image`.

Returns `(best_model, result::LMResult)`.

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

Pass `reweight` as a `LossFunctions.SupervisedLoss` (e.g. `LossFunctions.HuberLoss()`,
`PSFModels.TukeyLoss()`) to enable iteratively reweighted least squares (IRLS). After
each accepted LM step, the residuals are used to recompute pixel weights via
``w_i^{\\text{final}} = w_i^{\\text{base}} \\cdot w(r_i/\\sigma)`` where
``w(r) = \\psi(r)/r`` is derived from the loss function's influence function
and ``\\sigma`` is a robust scale estimate. This provides automatic outlier
rejection — useful for cosmic rays, bad pixels, and satellite trails in
astronomical images.

The scale ``\\sigma`` is estimated by `scale_estimator`. If not provided,
it defaults to `FixedScale` inferred from `inv_var` if supplied,
otherwise `MADScale()` is used. Available scale estimators:
- `MADScale()` — median absolute deviation (robust, default without `inv_var`)
- `FixedScale(σ)` — fixed user-provided scale
- `MScale(; δ=0.5)` — iterative M-scale (most robust, highest cost)

Recommended loss functions for IRLS:
- `LossFunctions.HuberLoss(c)` — soft downweighting of outliers, suitable for mildly
  contaminated data or crowded-field photometry
- `PSFModels.TukeyLoss(; c=4.685)` — complete rejection of extreme outliers, ideal
  for cosmic rays and bad pixels in space-based imaging

# Covariance Estimation
The covariance of the fitted parameters is estimated from the Gauss-Newton Hessian
approximation to the Hessian at the solution. For known good input weights `inv_var`,
the covariance is simply the inverse of this Hessian approximation. Use 
`covariance_estimator = KnownWeightsCovarianceEstimator()` for this behavior.
However, when IRLS reweighting is used and the weights are estimated from the data,
the covariance is inflated by the reduced cost per degree of freedom to account for
uncertainty in the weights. In this case, use `ReweightedCovarianceEstimator`.

# Damping Strategies
- `Marquardt`: diagonal entries are scaled by `max(A[i, i], min_diagonal)` to
  prevent small curvature directions from being under-damped
- `Levenberg`: uniform damping with no scaling
- `NoDamping`: no damping; equivalent to Gauss-Newton (not recommended)

# Keyword arguments

- `fixed`: `NamedTuple` of frozen parameter name → value pairs
- `inv_var`: inverse-variance weights, same shape as `image`
- `reweight`: a `LossFunctions.SupervisedLoss` for IRLS reweighting,
  or `nothing` (default) for standard L2 fitting
- `scale_estimator::AbstractScaleEstimator`: scale estimator for IRLS;
  defaults to `FixedScale` (from `inv_var`) or `MADScale()` (see above)
- `covariance_estimator::AbstractCovarianceEstimator`: method for estimating the covariance matrix of the fitted parameters; defaults to `ReweightedCovarianceEstimator` when IRLS is used and `KnownWeightsCovarianceEstimator` otherwise
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
"""
function fit_lm(
        model::AbstractPSFModel{T},
        image::AbstractMatrix,
        inds = axes(image);
        fixed::NamedTuple = (;),
        inv_var = nothing,
        λ_init::Real = 1.0e-4,
        λ_up::Real = 10.0,
        λ_down::Real = 10.0,
        λ_min::Real = 1.0e-12,
        λ_max::Real = 1.0e12,
        damping::AbstractLMDamping = MarquardtDamping(),
        max_iter::Integer = 200,
        x_tol::Real = 1.0e-8,
        f_tol::Real = 1.0e-8,
        g_tol::Real = 1.0e-8,
        show_trace::Bool = false,
        reweight::Union{Nothing, LossFunctions.SupervisedLoss} = nothing,
        scale_estimator::Union{Nothing, AbstractScaleEstimator} = nothing,
        covariance_estimator::Union{Nothing, AbstractCovarianceEstimator} = nothing
    ) where {T}

    # Validate inputs
    if !isnothing(inv_var)
        size(inv_var) == size(image) ||
            throw(ArgumentError("`inv_var` must be the same size as `image`"))
        all(x -> isfinite(x) && x > 0, inv_var) ||
            throw(ArgumentError("`inv_var` must be finite and > 0 everywhere"))

    end
    _has_deriv(model) ||
        throw(
        ArgumentError(
            "model does not implement `evaluate_fg`; " *
                "Levenberg-Marquardt requires gradient"
        )
    )

    # Converting here simplifies downstream indexing
    inds = CartesianIndices(inds)

    # Parameter bookkeeping
    free_names, free_idx, x0 = free_params(model, fixed)
    n = length(x0)
    n > 0 || throw(ArgumentError("all model parameters are fixed; nothing to fit"))
    dof = length(inds) - n
    dof > 0 || throw(
        ArgumentError(
            "degrees of freedom must be positive; " *
                "too many free parameters ($n) for the number of pixels ($(prod(length, inds)))"
        )
    )
    free_names_val = Val(free_names)

    # Pre-allocate for in-place accumulation
    FT = float(T)
    λ = FT(λ_init)
    residuals = zeros(FT, axes(image))
    x, x_cand = Vector{FT}(x0), Vector{FT}(undef, n)
    A, b = zeros(FT, n, n), zeros(FT, n)
    A_cand, b_cand = zeros(FT, n, n), zeros(FT, n)
    A_damp = zeros(FT, n, n)
    δ = zeros(FT, n)

    # IRLS setup
    do_irls = !isnothing(reweight)
    # Determine scale estimator: explicit > FixedScale from inv_var > MAD
    scale_est = if !isnothing(scale_estimator)
        scale_estimator
    elseif !isnothing(inv_var)
        # FixedScale from inv_var — σ ≈ 1/√(inv_var) averaged over the cutout
        FixedScale(sqrt(inv(mean(view(inv_var, inds)))))
    else
        MADScale()
    end
    σ_final = FT(NaN) # Scale estimate from IRLS (NaN if not applicable)
    # Pre-allocate weight arrays for IRLS
    weights = nothing # Declare so it can be updated in-place if do_irls
    if do_irls
        weights = similar(image, FT)
        isnothing(inv_var) ? fill!(weights, one(FT)) : copyto!(weights, inv_var)
    end

    # Initial cost and Jacobian
    m0 = model_from_vector(model, free_names_val, x, fixed)
    cost = _lm_accum_fg!(A, b, m0, free_idx, image, inds, residuals, inv_var, weights)
    cost_init = cost

    if show_trace
        gnorm0 = sqrt(sum(abs2, b))
        println("Initialization | cost = $cost_init | λ = $λ | ||g|| = $gnorm0")
    end

    # LM iteration
    x_converged = false
    f_converged = false
    g_converged = false
    converged = false
    iter = 0

    while iter < max_iter
        iter += 1

        # Test for gradient-norm convergence before attempting the linear solve,
        # since if the gradient is small we are already close to a local minimum
        # and the step direction may be unreliable
        gnorm = sqrt(sum(abs2, b))
        if gnorm ≤ g_tol
            g_converged = true
            converged = true
            show_trace && println("Iter $(lpad(iter, 4)) | converged on gradient norm (||g|| = $gnorm)")
            break
        end

        # Build the damped normal matrix A + λD
        A_damp .= A
        damp!(A_damp, damping, λ)

        # Solve (A + λD) δ = −b; same as δ = A_damp \ (-b)
        local F
        try
            F = cholesky!(Symmetric(A_damp))
        catch e
            # If the matrix is not positive definite, increase λ and retry
            e isa LinearAlgebra.PosDefException || rethrow()
            λ = min(λ * FT(λ_up), FT(λ_max))
            continue
        end
        ldiv!(δ, F, -b)
        δnorm = sqrt(sum(abs2, δ))

        # Evaluate the candidate step
        x_cand .= x .+ δ
        m_cand = model_from_vector(model, free_names_val, x_cand, fixed)
        cost_cand = _lm_accum_fg!(A_cand, b_cand, m_cand, free_idx, image, inds, residuals, inv_var, weights)
        accepted = cost_cand < cost

        if show_trace
            status = accepted ? "accepted" : "rejected"
            println(
                "Iter $(lpad(iter, 4)) | cost = $cost → $cost_cand | " *
                    "λ = $λ | ||g|| = $gnorm | ||δ|| = $δnorm | $status"
            )
        end

        if accepted
            Δcost = cost - cost_cand
            cost = cost_cand
            λ = max(λ / FT(λ_down), FT(λ_min))
            x .= x_cand
            A .= A_cand
            b .= b_cand

            # Check standard LM convergence before reweighting
            x_converged = δnorm ≤ x_tol * (sqrt(sum(abs2, x)) + x_tol)
            f_converged = Δcost ≤ f_tol * (abs(cost) + f_tol)

            if x_converged || f_converged
                converged = true
                break
            end

            # IRLS: recompute weights from current residuals
            if do_irls
                σ_final = estimate_scale(scale_est, view(residuals, inds))

                # Skip reweighting if scale is zero (noiseless or degenerate)
                if σ_final > eps(FT)
                    for idx in CartesianIndices(inds)
                        r_scaled = residuals[idx] / σ_final
                        ivar = isnothing(inv_var) ? one(FT) : inv_var[idx]
                        weights[idx] = ivar * weight(reweight, FT(r_scaled))
                    end
                    # Recompute cost and gradient with new weights
                    cost = _lm_accum_fg!(A, b, m_cand, free_idx, image, inds, residuals, inv_var, weights)
                    # Reset λ — the Gauss-Newton approximation changed with the new weights
                    λ = λ_init
                end
            end
        else
            λ = min(λ * FT(λ_up), FT(λ_max))
        end
    end

    best_model = model_from_vector(model, free_names_val, x, fixed)

    covariance_estimator = if !isnothing(covariance_estimator)
        covariance_estimator
    elseif !isnothing(inv_var) && !do_irls
        KnownWeightsCovarianceEstimator()
    else
        ReweightedCovarianceEstimator()
    end

    cov = covariance!(covariance_estimator, A, cost, dof)
    σ² = cost / dof # mean weighted squares residual per dof in data units
    χ² = if isnothing(inv_var) && do_irls && !isnan(σ_final) && σ_final > eps(FT)
        σ² / σ_final^2 # reduced chi-squared: cost per dof normalized by the estimated scale²
    else
        σ² # reduced chi-squared without scale normalization (e.g. when inv_var is provided or no IRLS)
    end

    result = LMResult(
        x, cost, cost_init,
        converged, x_converged, f_converged, g_converged,
        iter, λ, σ_final, cov, χ²
    )
    return best_model, result
end

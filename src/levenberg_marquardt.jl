
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
"""
struct LMResult{T,V<:AbstractVector{T}}
    minimizer::V
    minimum::T
    cost_init::T
    converged::Bool
    x_converged::Bool
    f_converged::Bool
    g_converged::Bool
    iterations::Int
    λ_final::T
end

function Base.show(io::IO, r::LMResult)
    println(io, "LMResult:")
    println(io, "  converged:  ", r.converged,
                "  (x: ", r.x_converged,
                ", f: ", r.f_converged,
                ", g: ", r.g_converged, ")")
    println(io, "  iterations: ", r.iterations)
    println(io, "  cost:       ", r.minimum, "  (init: ", r.cost_init, ")")
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
# fit_lm
# ---------------------------------------------------------------------------

# Build the Gauss-Newton approximation to the Hessian (A = Jᵀ W J) and the
# gradient (b = Jᵀ W r) together with the cost C = ∑ wᵢ rᵢ².
#
# The Jacobian row for pixel i is `g_full[free_idx]`, the projection of the
# full analytic gradient onto the free parameters.  g_full is a StaticArrays
# SVector (immutable), so we index into it rather than mutating it.
function _lm_accum_fg!(A::AbstractMatrix{FT}, b::AbstractVector{FT}, m::AbstractPSFModel, free_idx, image, inds, inv_var) where {FT}
    @assert size(A, 1) == size(A, 2) == length(b) == length(free_idx)
    fill!(A, zero(FT))
    fill!(b, zero(FT))
    cost = zero(FT)
    n = length(free_idx)
    for idx in CartesianIndices(inds)
        px, py = idx[1], idx[2]
        w = isnothing(inv_var) ? one(FT) : FT(inv_var[idx])
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
                     show_trace=false)

Fit the free parameters of `model` to `image[inds]` under weighted L2 loss
using the Levenberg-Marquardt algorithm.  The model must implement
`evaluate_fg`.

`fixed` is a `NamedTuple` of field-name → value pairs whose parameters are
frozen during the fit.  All other fields of `model` are free.

Inverse variance weights can be passed via `inv_var`; it must be the same size
as `image`.  When provided, the covariance matrix of the free parameters is
also returned.

Returns `(best_model, result::LMResult)`, or `(best_model, result, cov)` when
`inv_var` is provided.

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

# Keyword arguments

- `fixed`: `NamedTuple` of frozen parameter name → value pairs
- `inv_var`: inverse-variance weights, same shape as `image`
- `reweight`: if `true`, allows `inv_var` to be updated based on residuals (default `false`)
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
                reweight::Bool=false) where {T}

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

    # Initial cost and Jacobian
    m0 = model_from_vector(model, free_names_val, x, fixed)
    cost = _lm_accum_fg!(A, b, m0, free_idx, image, inds, inv_var)
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
        F = cholesky!(F_buffer)
        ldiv!(δ, F, -b)
        δnorm = sqrt(sum(abs2, δ))

        # Evaluate the candidate step
        x_cand .= x .+ δ
        m_cand = model_from_vector(model, free_names_val, x_cand, fixed)
        cost_cand = _lm_accum_fg!(A_cand, b_cand, m_cand, free_idx, image, inds, inv_var)
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

            x_converged = δnorm ≤ FT(x_tol) * (sqrt(sum(abs2, x)) + FT(x_tol))
            f_converged = Δcost ≤ FT(f_tol) * (abs(cost) + FT(f_tol))

            if x_converged || f_converged
                converged = true
                break
            end
        else
            λ = min(λ * FT(λ_up), FT(λ_max))
        end
    end

    best_model = model_from_vector(model, free_names_val, x, fixed)
    result = LMResult(x, cost, cost_init,
                      converged, x_converged, f_converged, g_converged,
                      iter, λ)

    """
        compute_cov!(JTJ, cost, n_obs, n_params, inv_var, reweight::Bool)
        
    Compute the covariance matrix of the free parameters at the solution from the Gauss-Newton approximation to the Hessian (JTJ = Jᵀ W J for weights W) and the cost (which must be the weighted sum of squares, Χ² = rᵀ W r) at the solution. If `!isnothing(inv_var)`, then the covariance is given by the Cramér-Rao bound `inv(Jᵀ W J)`. Otherwise, we estimate the noise variance `σ² = cost / (n_obs - n_params)` and return `σ² * inv(Jᵀ J)`, which is the standard covariance estimate for nonlinear least squares and is asymptotically unbiased when the model is correct. If the relative weighting encoded by `inv_var` is trusted but its overall normalization may be inaccurate, you can set `reweight=true` to re-scale the covariance by the factor `cost / dof` to account for any mis-estimation of the noise variance.
    """
    function compute_cov!(JTJ, cost, n_obs, n_params, inv_var, reweight::Bool)
        F = cholesky!(JTJ) # modfies JTJ in-place to compute the Cholesky factorization
        cov = F \ I # = inv(JTJ), more stable
        if isnothing(inv_var) || reweight
            dof = n_obs - n_params
            dof > 0 || throw(ArgumentError("degrees of freedom must be positive to compute covariance"))
            σ² = cost / dof
            cov *= σ²
        end
        return cov
    end

    cov = compute_cov!(A, cost, prod(length, inds), n, inv_var, reweight)
    return best_model, result, cov
end

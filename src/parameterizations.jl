"""
    sum_massfractions(m)

Sum the mass-fraction weights across all supplied `MassFraction`s on their grid.
Returns a `Field` whose `tracer` is the per-cell total; dry points are `NaN`.
"""
function sum_massfractions(m::Union{NamedTuple,Vector})
    firstm = first(m)
    γ = firstm.γ
    T = eltype(firstm.fraction)
    nint = sum(γ.interior)
    ngrid = size(γ.wet)

    total = zeros(T, ngrid)
    Iint = cartesianindex(γ.interior)
    for m1 in m
        sum(m1.γ.interior) == nint || error("mass fraction grids do not match")
    end
    for I in Iint
        s = zero(T)
        for m1 in m
            m1.γ.wet[I] || continue
            s += m1.fraction[I]
        end
        total[I] = s
    end
    total[.!γ.wet] .= zero(T)/zero(T) # NaN on dry points

    return Field(total,
        γ,
        :sum_massfraction,
        "sum of mass fractions",
        "unitless")
end

"""
    softmax_vector(v; α = 1)

Compute a (temperatured) softmax of vector `v`, producing non-negative components
that sum to 1. Uses log-sum-exp for numerical stability. Larger `α` flattens; `α=1`
is the standard softmax.
"""
function softmax_vector(v::Vector{T}; α::Real = 1) where T
    isempty(v) && return T[]
    vmax = maximum(v ./ α)
    shifted = (v ./ α .- vmax)
    logden = log(sum(exp.(shifted)))
    s = exp.(shifted .- logden)
    return s
end

"""
    softmax_vector!(v; α = 1)

In-place softmax of `v`. Overwrites `v` with its temperatured softmax (non-negative,
sum to 1) using a stable log-sum-exp. Returns `v`. Works with vectors and vector
views.
"""
function softmax_vector!(v::AbstractVector{T}; α::Real = 1) where T
    isempty(v) && return v
    invα = inv(α)
    vmax = @inbounds v[1] * invα
    @inbounds for i in 2:length(v)
        val = v[i] * invα
        val > vmax && (vmax = val)
    end

    denom = zero(T)
    @inbounds for i in eachindex(v)
        v[i] = exp(v[i] * invα - vmax)
        denom += v[i]
    end
    invden = inv(denom)
    @inbounds for i in eachindex(v)
        v[i] *= invden
    end
    return v
end


"""
    gsoftmax_vector(gs, s; α = 1)

Adjoint of `softmax_vector` with temperature `α`. Gradients scale like `1/α`.
Works with vectors and vector views.
"""
function gsoftmax_vector(gs::AbstractVector{T}, s::AbstractVector{T}; α::Real = 1) where T
    work_gm = similar(s)
    gsoftmax_vector!(work_gm, gs, s; α = α)
    return work_gm
end

"""
    gsoftmax_vector!(gv, gs, s; α = 1)

In-place adjoint of `softmax_vector`. Writes gradients into `gv` (matching length),
given upstream `gs` and softmaxed `s`. Works with vectors and views. Returns `gv`.
"""
function gsoftmax_vector!(work_gm::AbstractVector{T}, gs::AbstractVector{T}, s::AbstractVector{T}; α::Real = 1) where T
    (length(s) == length(gs) && length(gs) == length(work_gm)) || error("gsoftmax_vector!: length mismatch")
    dotgs = dot(gs, s)
    scale = inv(α)
    @inbounds for i in eachindex(s)
        work_gm[i] = scale * s[i] * (gs[i] - dotgs)
    end
    return work_gm
end

"""
    softmax_massfractions(m; α = 1)

Apply a standard softmax across neighboring
mass fractions at each grid cell and
return a new vector of `MassFraction`s with fractions that sum to 1 per cell.
Dry points remain `NaN`.
"""
function softmax_massfractions(m::NamedTuple; α::Real = 1)
    firstm = first(m)
    γ = firstm.γ
    T = eltype(firstm.fraction)
    nint = sum(γ.interior)
    ngrid = size(γ.wet)

    for k in eachindex(m)
        sum(m[k].γ.interior) == nint || error("mass fraction grids do not match")
    end

    n_neighbors = neighbors(m, γ)
    max_neighbors = maximum(n_neighbors.tracer)
    work = Vector{T}(undef, max(max_neighbors, 1))

    transformed_m = similar(m)

    for I in eachindex(n_neighbors.tracer)
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue
        j = 1
        for k in eachindex(m)
            if m[k].γ.wet[I]
                work[j] = m[k].fraction[I]
                j += 1
            end
        end

        softmax_vector!(view(work, 1:nneigh); α = α)

        j = 1
        for k in eachindex(m)
            if m[k].γ.wet[I]
                transformed_m[k].fraction[I] = work[j]
                j += 1
            end
        end
    end
    s = transformed_m
    return s
end

"""
    softmax_massfractions!(x; α = 1)

In-place version of `softmax_massfractions`. Softmaxes `x` on each wet cell and
writes the results back into `x` without allocating new `MassFraction`s. Returns `x`.
"""
function softmax_massfractions!(x::NamedTuple; α::Real = 1)
    firstm = first(x)
    γ = firstm.γ
    T = eltype(firstm.fraction)
    nint = sum(γ.interior)

    @inbounds for k in eachindex(x)
        sum(x[k].γ.interior) == nint || error("softmax_massfractions!: grid mismatch")
    end

    n_neighbors = neighbors(x, γ)
    max_neighbors = maximum(n_neighbors.tracer)
    work = Vector{T}(undef, max(max_neighbors, 1))

    for I in eachindex(n_neighbors.tracer)
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue

        j = 1
        for k in eachindex(x)
            if x[k].γ.wet[I]
                work[j] = x[k].fraction[I]
                j += 1
            end
        end

        softmax_vector!(view(work, 1:nneigh); α = α)

        j = 1
        for k in eachindex(x)
            if x[k].γ.wet[I]
                x[k].fraction[I] = work[j]
                j += 1
            end
        end
    end

    # Ensure dry points stay NaN to keep vec() consistent with wet mask.
    nanval = zero(T) / zero(T)
    # @inbounds for k in eachindex(x)
    #     wet = x[k].γ.wet
    #     frac = x[k].fraction
    #     for I in eachindex(wet)
    #         wet[I] && continue
    #         frac[I] = nanval
    #     end
    # end

    return x
end


"""
    gsoftmax_massfractions(gs, s; α = 1)

Adjoint of `softmax_massfractions`: given softmaxed mass fractions `s` and upstream
sensitivities `gs` (matching NamedTuple layout), return gradients with respect to the
pre-softmax mass fractions. Uses `gsoftmax_vector` on each interior cell.
"""
function gsoftmax_massfractions(gs::NamedTuple, s::NamedTuple; α::Real = 1)
    keys(s) == keys(gs) || error("gsoftmax_massfractions: key mismatch")
    γ = first(s).γ
    T = eltype(first(s).fraction)
    nint = sum(γ.interior)

    for k in keys(s)
        sum(s[k].γ.interior) == nint || error("gsoftmax_massfractions: grid mismatch")
    end

    # initialize gradient container; fill everything with NaN to make
    # unintended use obvious and keep non-wet points out of vec().
    g_m = map(s) do mf
        gmf = similar(mf)
        gmf.fraction .= zero(T) / zero(T) # NaN everywhere
        gmf
    end

    n_neighbors = neighbors(s, γ)
    max_neighbors = maximum(n_neighbors.tracer)
    work_s = Vector{T}(undef, max(max_neighbors, 1))
    work_gs = Vector{T}(undef, max(max_neighbors, 1))
    work_gm = Vector{T}(undef, max(max_neighbors, 1))

    for I in eachindex(n_neighbors.tracer)
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue
        j = 1
        for k in eachindex(s)
            if s[k].γ.wet[I]
                work_s[j] = s[k].fraction[I]
                work_gs[j] = gs[k].fraction[I]
                j += 1
            end
        end
        gvec = gsoftmax_vector!(view(work_gm, 1:nneigh),
                                view(work_gs, 1:nneigh),
                                view(work_s, 1:nneigh); α = α)

        j = 1
        for k in eachindex(s)
            if s[k].γ.wet[I]
                g_m[k].fraction[I] = gvec[j]
                j += 1
            end
        end
    end

    return g_m
end

"""
    gsoftmax_massfractions!(gs, s; α = 1)

In-place version of `gsoftmax_massfractions`. Overwrites `gs` with gradients of the
pre-softmax mass fractions given upstream sensitivities `gs` and softmaxed mass
fractions `s`. Returns `gs`.
"""
function gsoftmax_massfractions!(gs::NamedTuple, s::NamedTuple; α::Real = 1)
    keys(s) == keys(gs) || error("gsoftmax_massfractions!: key mismatch")
    γ = first(s).γ
    T = eltype(first(s).fraction)
    nint = sum(γ.interior)

    for k in keys(s)
        sum(s[k].γ.interior) == nint || error("gsoftmax_massfractions!: grid mismatch")
        sum(gs[k].γ.interior) == nint || error("gsoftmax_massfractions!: grid mismatch")
        gs[k].fraction .= zero(T) / zero(T) # set to NaN everywhere
    end

    n_neighbors = neighbors(s, γ)
    max_neighbors = maximum(n_neighbors.tracer)
    work_s = Vector{T}(undef, max(max_neighbors, 1))
    work_gs = Vector{T}(undef, max(max_neighbors, 1))
    work_gm = Vector{T}(undef, max(max_neighbors, 1))

    for I in eachindex(n_neighbors.tracer)
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue
        j = 1
        for k in eachindex(s)
            if s[k].γ.wet[I]
                work_s[j] = s[k].fraction[I]
                work_gs[j] = gs[k].fraction[I]
                j += 1
            end
        end

        gvec = gsoftmax_vector!(view(work_gm, 1:nneigh),
                                view(work_gs, 1:nneigh),
                                view(work_s, 1:nneigh); α = α)

        j = 1
        for k in eachindex(s)
            if s[k].γ.wet[I]
                gs[k].fraction[I] = gvec[j]
                j += 1
            end
        end
    end

    gm = gs #gs has now been transformed to gm
    
    return gm
end


"""
    invsoftmax_vector(s; α = 1)

Convert softmax outputs back to a pre-softmax vector. Softmax ignores any constant
offset, so we return the simple choice `log.(s)`, which re-softmaxes to `s` when
entries are positive and sum to 1.
"""
function invsoftmax_vector(s::Vector{T}; α::Real = 1) where T
    any(x -> !(x > zero(T)), s) && error("invsoftmax_vector: inputs must be positive")
    return α .* log.(s)
end


"""
    invsoftmax_massfractions(s; α = 1)

Inverse of `softmax_massfractions`. Takes softmaxed mass fractions `s` and logs each
wet point, producing pre-softmax values that softmax back to `s`. Dry points stay `NaN`.
"""
function invsoftmax_massfractions(s::NamedTuple; α::Real = 1)
    γ = first(s).γ
    T = eltype(first(s).fraction)
    nint = sum(γ.interior)

    for k in keys(s)
        sum(s[k].γ.interior) == nint || error("invsoftmax_massfractions: grid mismatch")
    end

    log_m = similar(s)
    n_neighbors = neighbors(s, γ)

    Iint = cartesianindex(γ.interior)
    for I in Iint
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue
        for k in keys(s)
            if s[k].γ.wet[I]
                val = s[k].fraction[I]
                val > zero(T) || error("invsoftmax_massfractions: softmax components must be positive")
                log_m[k].fraction[I] = α * log(val)
            end
        end
    end

    return log_m
end

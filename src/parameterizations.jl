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
    gsoftmax_vector(gs, s; α = 1)

Adjoint of `softmax_vector` with temperature `α`. Gradients scale like `1/α`.
"""
function gsoftmax_vector(gs::Vector{T}, s::Vector{T}; α::Real = 1) where T
    length(s) == length(gs) || error("gs and v have length mismatch")
    gv = (1/α) * s .* (gs .- dot(gs, s))
    return gv
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

    m_keys = collect(keys(m))
    nf = length(m_keys)

    for k in m_keys
        sum(m[k].γ.interior) == nint || error("mass fraction grids do not match")
    end

    n_neighbors = neighbors(m, γ)

    transformed_m = similar(m)

    Iint = cartesianindex(γ.interior)
    for I in Iint
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue
        untransformed_fractions = Vector{T}(undef, nneigh)
        j = 1
        for k in m_keys
            if m[k].γ.wet[I]
                untransformed_fractions[j] = m[k].fraction[I]
                j += 1
            end
        end
        transformed_fractions = softmax_vector(untransformed_fractions; α = α)

        j = 1
        for k in m_keys
            if m[k].γ.wet[I]
                transformed_m[k].fraction[I] = transformed_fractions[j]
                j += 1
            end
        end
    end
    s = transformed_m
    return s
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

    Iint = cartesianindex(γ.interior)
    for I in Iint
        nneigh = n_neighbors.tracer[I]
        nneigh == 0 && continue
        svec = Vector{T}(undef, nneigh)
        gsvec = Vector{T}(undef, nneigh)
        j = 1
        for k in keys(s)
            if s[k].γ.wet[I]
                svec[j] = s[k].fraction[I]
                gsvec[j] = gs[k].fraction[I]
                j += 1
            end
        end
        gvec = gsoftmax_vector(gsvec, svec; α = α)

        j = 1
        for k in keys(s)
            if s[k].γ.wet[I]
                g_m[k].fraction[I] = gvec[j]
                j += 1
            end
        end
    end

    return g_m
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

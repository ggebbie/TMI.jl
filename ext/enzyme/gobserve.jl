"""
    augmented_primal(config, observe, Duplicated, c, locs, γ)

Forward pass for point observations `y = observe(c, locs, γ)`.

The primal observation is computed from interpolation weights directly. The
reverse pass reuses the taped weights and applies the same adjoint pattern as
`gobserve` in `TMI.jl`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(observe)},
    ::Type{<:Duplicated},
    c::Duplicated{<:Field{T}},
    locs::Const{<:Vector{<:Tuple}},
    γ::Const{<:Grid},
) where {T <: Real}
    y = Vector{T}(undef, length(locs.val))
    wetmask = c.val.γ.wet
    weights = Vector{typeof(interpweights(first(locs.val), c.val.γ))}(undef, length(locs.val))
    for ii in eachindex(y)
        weights[ii] = interpweights(locs.val[ii], c.val.γ)
        y[ii] = sum(c.val.tracer[wetmask] .* weights[ii][wetmask])
    end
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    return AugmentedReturn(primal, gy, (gy, weights))
end

"""
    reverse(config, observe, Duplicated, tape, c, locs, γ)

Reverse pass for point observations.

Here `gy = dJ/dy`. Each observation contributes its interpolation weights back
to the tracer field, so the rule accumulates `gc = dJ/dc` into `c.dval`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(observe)},
    ::Type{<:Duplicated},
    tape,
    c::Duplicated{<:Field},
    locs::Const{<:Vector{<:Tuple}},
    γ::Const{<:Grid},
)
    gy, weights = tape
    for ii in eachindex(gy)
        c.dval.tracer[c.val.γ.wet] .+= gy[ii] * weights[ii][c.val.γ.wet]
    end

    return (nothing, nothing, nothing)
end

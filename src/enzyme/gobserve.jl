"""
    augmented_primal(config, observe, Duplicated, c, locs, γ)

Forward pass for point observations `y = observe(c, locs, γ)`.

The primal observation uses TMI's interpolation machinery. For the reverse pass
we tape the interpolation weights computed by `interpweights`, so `reverse` can
apply the same adjoint pattern as `gobserve` in `TMI.jl`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(observe)},
    ::Type{<:Duplicated},
    c::Annotation{<:Field},
    locs::Annotation{<:Vector{<:Tuple}},
    γ::Annotation{<:Grid},
)
    y = func.val(c.val, locs.val, γ.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    weights = needs_shadow(config) ? [interpweights(loc, c.val.γ) for loc in locs.val] : nothing

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
    c::Annotation{<:Field},
    locs::Annotation{<:Vector{<:Tuple}},
    γ::Annotation{<:Grid},
)
    gy, weights = tape
    if c isa Duplicated
        for ii in eachindex(gy)
            c.dval.tracer[c.val.γ.wet] .+= gy[ii] * weights[ii][c.val.γ.wet]
        end
    end

    return (nothing, nothing, nothing)
end

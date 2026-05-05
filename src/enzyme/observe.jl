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
    g_y = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    weights = needs_shadow(config) ? [interpweights(loc, c.val.γ) for loc in locs.val] : nothing

    return AugmentedReturn(primal, g_y, (g_y, weights))
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(observe)},
    ::Type{<:Duplicated},
    tape,
    c::Annotation{<:Field},
    locs::Annotation{<:Vector{<:Tuple}},
    γ::Annotation{<:Grid},
)
    g_y, weights = tape
    if c isa Duplicated
        for ii in eachindex(g_y)
            c.dval.tracer[c.val.γ.wet] .+= g_y[ii] * weights[ii][c.val.γ.wet]
        end
    end

    return (nothing, nothing, nothing)
end

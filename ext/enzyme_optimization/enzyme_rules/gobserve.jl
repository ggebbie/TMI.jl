"""
    augmented_primal(config, observe, return_activity, field,
        interpolation_indices, γ)
    reverse(config, observe, return_activity, tape, field,
        interpolation_indices, γ)

Evaluate sparse observations and reverse their interpolation.

# Arguments
- `config`, `field`: reverse configuration and duplicated field
- `interpolation_indices`, `γ`, `tape`: interpolation metadata and cotangent tape

# Output
- `result`: augmented samples or reverse placeholders
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(observe)},
    ::Type{<:Union{Duplicated,Enzyme.DuplicatedNoNeed}},
    field::Duplicated{<:Field{T}},
    interpolation_indices::Annotation{<:Vector{<:Tuple}},
    γ::Annotation{<:Grid},
) where {T<:Real}
    y = func.val(field.val, interpolation_indices.val, γ.val)
    primal = needs_primal(config) ? y : nothing
    g_y = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    wet_weights = needs_shadow(config) ?
        [interpweights(indices, γ.val)[field.val.γ.wet]
            for indices in interpolation_indices.val] : nothing
    return AugmentedReturn(primal, g_y, (g_y, wet_weights))
end
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(observe)},
    ::Type{<:Union{Duplicated,Enzyme.DuplicatedNoNeed}},
    tape,
    field::Duplicated{<:Field},
    ::Annotation{<:Vector{<:Tuple}},
    ::Annotation{<:Grid},
)
    g_y, wet_weights = tape
    wet_mask = field.val.γ.wet
    for sample_index in eachindex(g_y)
        field.dval.tracer[wet_mask] .+=
            g_y[sample_index] .* wet_weights[sample_index]
    end
    return (nothing, nothing, nothing)
end

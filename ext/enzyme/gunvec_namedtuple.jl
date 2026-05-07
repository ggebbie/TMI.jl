"""
    augmented_primal(config, unvec, Duplicated, template, uvec)
    augmented_primal(config, unvec, MixedDuplicated, template, uvec)

Forward rules for

    y = unvec(template::NamedTuple, uvec)

where `template` is a `NamedTuple` whose values are tracer-like TMI structs
(`Source`, `Field`, `BoundaryCondition`, `MassFraction`). The rule returns `y`
and allocates `gy = dJ/dy` for the reverse pass.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{NT},
    uvec::Duplicated{Vector{T}},
) where {names, U <: Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}, NT <: NamedTuple{names, U}, T <: Real}
    y = needs_primal(config) ? func.val(template.val, uvec.val) : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:MixedDuplicated},
    template::Const{NT},
    uvec::Duplicated{Vector{T}},
) where {names, U <: Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}, NT <: NamedTuple{names, U}, T <: Real}
    y = func.val(template.val, uvec.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    shadow = isnothing(gy) ? nothing : Ref(gy)
    return AugmentedReturn(primal, shadow, shadow)
end

"""
    reverse(config, unvec, Duplicated, gy, template, uvec)
    reverse(config, unvec, MixedDuplicated, gy, template, uvec)

Reverse rules for `unvec(::NamedTuple, ::Vector)`.

Each tuple entry contributes a contiguous segment of `uvec`; the reverse rule
accumulates entry-wise wet gradients from `gy` back into `uvec.dval` using the
same packed layout as `vec(template)`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    gy::NT,
    template::Const{NT},
    uvec::Duplicated{Vector{T}},
) where {names, U <: Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}, NT <: NamedTuple{names, U}, T <: Real}
    nlo = 1
    for (gyi, ti) in zip(gy, template.val)
        n = sum(wet(ti))
        uvec.dval[nlo:nlo+n-1] .+= gyi isa MassFraction ? gyi.fraction[wet(ti)] : gyi.tracer[wet(ti)]
        nlo += n
    end
    return (nothing, nothing)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:MixedDuplicated},
    gy::Base.RefValue{NT},
    template::Const{NT},
    uvec::Duplicated{Vector{T}},
) where {names, U <: Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}, NT <: NamedTuple{names, U}, T <: Real}
    nlo = 1
    for (gyi, ti) in zip(gy[], template.val)
        n = sum(wet(ti))
        uvec.dval[nlo:nlo+n-1] .+= gyi isa MassFraction ? gyi.fraction[wet(ti)] : gyi.tracer[wet(ti)]
        nlo += n
    end
    return (nothing, nothing)
end

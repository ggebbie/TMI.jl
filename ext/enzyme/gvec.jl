function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(vec)},
    ::Type{<:Duplicated},
    u::Annotation{NT},
) where {
    names,
    NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}},
}
    y = needs_primal(config) ? func.val(u.val) : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(u.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(vec)},
    ::Type{<:MixedDuplicated},
    u::Annotation{NT},
) where {
    names,
    NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}},
}
    y = func.val(u.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    shadow = isnothing(gy) ? nothing : Ref(gy)
    return AugmentedReturn(primal, shadow, shadow)
end

"""
    reverse(config, vec, Duplicated/MixedDuplicated, gy, u)

Adjoint for `vec`: scatter vector gradients back into tracer storage.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(vec)},
    ::Type{<:Duplicated},
    gy::Vector{T},
    u::Annotation{NT},
) where {
    names,
    T<:Real,
    NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}},
}
    if !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        nlo = 1
        for (udi, ui) in zip(ud, u.val)
            n = sum(wet(ui))
            udi.tracer[wet(ui)] .+= gy[nlo:nlo+n-1]
            nlo += n
        end
    end
    return (nothing,)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(vec)},
    ::Type{<:MixedDuplicated},
    gy::Base.RefValue{Vector{T}},
    u::Annotation{NT},
) where {
    names,
    T<:Real,
    NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}},
}
    if !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        nlo = 1
        for (udi, ui) in zip(ud, u.val)
            n = sum(wet(ui))
            udi.tracer[wet(ui)] .+= gy[][nlo:nlo+n-1]
            nlo += n
        end
    end
    return (nothing,)
end

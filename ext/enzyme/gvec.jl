"""
    augmented_primal(config, vec, Duplicated/MixedDuplicated, u)

Custom `vec` rules for heterogeneous NamedTuple controls. This avoids Enzyme
type-activity ambiguity on packed controls used by the inversion test path.

Linked TMI function:
- `TMI.vec`

Arguments:
- `config`: Enzyme reverse-mode config.
- `func`: wrapped `vec`.
- `u`: control NamedTuple annotation.
"""
function augmented_primal(config::RevConfigWidth{1}, func::Const{typeof(vec)},
    ::Type{<:Duplicated}, u::Annotation{NT}) where
    {names, NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}}}
    y = needs_primal(config) ? func.val(u.val) : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(u.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

function augmented_primal(config::RevConfigWidth{1}, func::Const{typeof(vec)},
    ::Type{<:MixedDuplicated}, u::Annotation{NT}) where
    {names, NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}}}
    y = func.val(u.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    shadow = isnothing(gy) ? nothing : Ref(gy)
    return AugmentedReturn(primal, shadow, shadow)
end

"""
    reverse(config, vec, Duplicated/MixedDuplicated, gy, u)

Adjoint for `vec`: scatter vector gradients back into control storage
(`.tracer` for tracer-like controls, `.fraction` for MassFraction).

Linked TMI function:
- `TMI.vec`

Arguments:
- `gy`: gradient w.r.t. packed vector output.
- `u`: control NamedTuple annotation with shadow storage.
"""
function reverse(::RevConfigWidth{1}, ::Const{typeof(vec)}, ::Type{<:Duplicated},
    gy::Vector{T}, u::Annotation{NT}) where
    {names, T<:Real, NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}}}
    if !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        nlo = 1
        for (udi, ui) in zip(ud, u.val)
            n = sum(wet(ui))
            if udi isa MassFraction
                udi.fraction[wet(ui)] .+= gy[nlo:nlo+n-1]
            else
                udi.tracer[wet(ui)] .+= gy[nlo:nlo+n-1]
            end
            nlo += n
        end
    end
    return (nothing,)
end

function reverse(::RevConfigWidth{1}, ::Const{typeof(vec)}, ::Type{<:MixedDuplicated},
    gy::Base.RefValue{Vector{T}}, u::Annotation{NT}) where
    {names, T<:Real, NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition, MassFraction}}}}}
    if !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        nlo = 1
        for (udi, ui) in zip(ud, u.val)
            n = sum(wet(ui))
            if udi isa MassFraction
                udi.fraction[wet(ui)] .+= gy[][nlo:nlo+n-1]
            else
                udi.tracer[wet(ui)] .+= gy[][nlo:nlo+n-1]
            end
            nlo += n
        end
    end
    return (nothing,)
end

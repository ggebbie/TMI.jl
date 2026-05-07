"""
    augmented_primal(config, unvec!, Const, u, uvec)

Forward rule for in-place unpacking of compact controls:

    unvec!(u, uvec)

Supports tracer-like `u` and NamedTuples of tracer-like TMI structs.
"""
function augmented_primal(
    ::RevConfigWidth{1},
    func::Const{typeof(unvec!)},
    ::Type{<:Const},
    u::Annotation{U},
    uvec::Annotation{<:Vector{T}},
) where {U<:Union{Field, BoundaryCondition, Source}, T<:Real}
    func.val(u.val, uvec.val)
    # `unvec!` is mutating and returns `nothing` in Julia, so the primal return
    # in Enzyme's AugmentedReturn must also be `nothing` (not a rebuilt `y`).
    return AugmentedReturn(nothing, nothing, nothing)
end

function augmented_primal(
    ::RevConfigWidth{1},
    func::Const{typeof(unvec!)},
    ::Type{<:Const},
    u::Annotation{NT},
    uvec::Annotation{<:Vector{T}},
) where {
    names,
    NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}},
    T<:Real,
}
    func.val(u.val, uvec.val)
    # Same reasoning here: gradients flow through mutated `u`/`u.dval`; there is
    # no returned object to attach `gy` to for this in-place API.
    return AugmentedReturn(nothing, nothing, nothing)
end

"""
    reverse(config, unvec!, Const, tape, u, uvec)

Adjoint for `unvec!`: gather gradients from rebuilt controls back into `uvec`
with the same packing as `vec`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec!)},
    ::Type{<:Const},
    tape,
    u::Annotation{U},
    uvec::Annotation{<:Vector{T}},
) where {U<:Union{Field, BoundaryCondition, Source}, T<:Real}
    if !(uvec isa Const) && !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        uvec.dval .+= ud.tracer[wet(u.val)]
    end
    return (nothing, nothing)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec!)},
    ::Type{<:Const},
    tape,
    u::Annotation{NT},
    uvec::Annotation{<:Vector{T}},
) where {
    names,
    NT<:NamedTuple{names, <:Tuple{Vararg{<:Union{Source, Field, BoundaryCondition}}}},
    T<:Real,
}
    if !(uvec isa Const) && !(u isa Const)
        ud = u isa MixedDuplicated ? u.dval[] : u.dval
        nlo = 1
        for (udi, ui) in zip(ud, u.val)
            n = sum(wet(ui))
            uvec.dval[nlo:nlo+n-1] .+= udi.tracer[wet(ui)]
            nlo += n
        end
    end
    return (nothing, nothing)
end

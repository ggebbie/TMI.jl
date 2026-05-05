"""
    augmented_primal(config, unvec, Duplicated, template, uvec)
    augmented_primal(config, unvec, MixedDuplicated, template, uvec)

Forward rule for reverse-mode AD through

    y = unvec(template::BoundaryCondition, uvec)

`unvec` is the inverse layout operation to `vec`: it rebuilds a full
`BoundaryCondition` from the compact vector of wet boundary values, using the
wet mask and metadata stored in `template`.

A shadow is Enzyme's storage for an adjoint value. Here, `gy` denotes the output
shadow associated with `y`, i.e. `gy = dJ/dy`.

Two methods are provided because `unvec` returns a `BoundaryCondition`, but the
derivative of that return value is not always handed back to the reverse rule in
the same form. `Duplicated` gives the reverse rule the derivative object
directly, while `MixedDuplicated` gives it a reference to that object.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    y  = needs_primal(config) ? func.val(template.val, uvec.val) : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:MixedDuplicated},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    y = func.val(template.val, uvec.val)

    primal = needs_primal(config) ? y : nothing
    gy     = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    shadow = isnothing(gy) ? nothing : Ref(gy)

    return AugmentedReturn(primal, shadow, shadow)
end

"""
    reverse(config, unvec, Duplicated, gy, template, uvec)
    reverse(config, unvec, MixedDuplicated, gy, template, uvec)

Reverse rule for

    y = unvec(template::BoundaryCondition, uvec)

`unvec` rebuilds a full `BoundaryCondition` from the compact vector `uvec`,
placing those values at the wet boundary entries defined by `template`.

A shadow is Enzyme's storage for an adjoint value.
Here, `gy` is the derivative storage for `y`, and `uvec.dval` is
the derivative storage for `uvec`.

The reverse rule reads the wet tracer entries of `gy` and accumulates them into
`uvec.dval`; equivalently, if `guvec = dJ/duvec` and `gy = dJ/dy`, then
`guvec += vec(gy)`, where `vec(gy)` denotes `gy.tracer[template.wet]`.

For `Duplicated`, the reverse rule reads directly from `gy`. For
`MixedDuplicated`, it reads from `gy[]`.

The accumulation uses `.+=` because other reverse paths may already have
contributed to `uvec.dval`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    gy::BoundaryCondition{T,R,N,G,B},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    uvec.dval .+= gy.tracer[template.val.wet]
    return (nothing, nothing)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:MixedDuplicated},
    gy::Base.RefValue{BoundaryCondition{T,R,N,G,B}},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    uvec.dval .+= gy[].tracer[template.val.wet]
    return (nothing, nothing)
end
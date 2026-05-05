"""
    augmented_primal(config, unvec, Duplicated, template, uvec)
    augmented_primal(config, unvec, MixedDuplicated, template, uvec)

Forward rule for reverse-mode AD through

    y = unvec(template::BoundaryCondition, uvec)

`unvec` is the inverse layout operation to `vec`: it rebuilds a full
`BoundaryCondition` from the compact vector of wet boundary values, using the
wet mask and metadata stored in `template`.

In Enzyme terminology, a shadow is the storage paired with an active primal
value. In reverse mode, the output shadow stores the adjoint accumulated from
downstream computations. Here, the output shadow `gy` represents `ȳ = ∂J/∂y`.

Two methods are needed because Enzyme uses different shadow conventions for
`Duplicated` and `MixedDuplicated` outputs. For `Duplicated`, the output shadow
is passed directly to the reverse rule as a `BoundaryCondition`. For
`MixedDuplicated`, the output shadow is passed by reference, so the reverse rule
receives a `Ref{BoundaryCondition}`. The mathematical rule is identical in both
cases; only the shadow container differs.
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

Since `unvec` scatters the compact vector `uvec` into the wet entries of a full
`BoundaryCondition`, its adjoint gathers the wet entries of the output shadow
back into the input shadow `uvec.dval`.

Equivalently, the adjoint of `unvec` is the corresponding boundary-condition
packing operation:

    ūvec += vec(ȳ)

where `ȳ` is the output shadow. This is implemented by accumulating the wet
tracer entries of `gy` into `uvec.dval`.

The two methods differ only in how Enzyme passes the output shadow. For
`Duplicated`, `gy` is received directly as a `BoundaryCondition`. For
`MixedDuplicated`, `gy` is received as a `Ref{BoundaryCondition}` and must be
accessed with `gy[]`.

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
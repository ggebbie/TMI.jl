"""
    augmented_primal(config, unvec, Duplicated, template, uvec)
    augmented_primal(config, unvec, MixedDuplicated, template, uvec)

Forward rules for

    y = unvec(template, uvec)

where `template` is one of `Field`, `BoundaryCondition`, or `Source`.

`unvec` rebuilds a tracer-like object from the compact control vector `uvec`
using `wet(template)` to define active entries. The rule returns `y` and
allocates the output shadow `gy = dJ/dy` for the reverse pass.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{U},
    uvec::Duplicated{Vector{T}},
) where {U<:Union{Field, BoundaryCondition, Source}, T<:Real}
    y = needs_primal(config) ? func.val(template.val, uvec.val) : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:MixedDuplicated},
    template::Const{U},
    uvec::Duplicated{Vector{T}},
) where {U<:Union{Field, BoundaryCondition, Source}, T<:Real}
    y = func.val(template.val, uvec.val)
    primal = needs_primal(config) ? y : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(y) : nothing
    shadow = isnothing(gy) ? nothing : Ref(gy)
    return AugmentedReturn(primal, shadow, shadow)
end

"""
    reverse(config, unvec, Duplicated, gy, template, uvec)
    reverse(config, unvec, MixedDuplicated, gy, template, uvec)

Reverse rules for tracer-like `unvec`.

The adjoint of `unvec` is masked `vec`: gradients from the rebuilt object are
read at `wet(template)` and accumulated into `uvec.dval`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    gy::U,
    template::Const{U},
    uvec::Duplicated{Vector{T}},
) where {U<:Union{Field, BoundaryCondition, Source}, T<:Real}
    uvec.dval .+= gy.tracer[wet(template.val)]
    return (nothing, nothing)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:MixedDuplicated},
    gy::Base.RefValue{U},
    template::Const{U},
    uvec::Duplicated{Vector{T}},
) where {U<:Union{Field, BoundaryCondition, Source}, T<:Real}
    uvec.dval .+= gy[].tracer[wet(template.val)]
    return (nothing, nothing)
end

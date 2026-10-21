"""
    function augmented_primal(config, unvec!, return_activity, destination, x)
    function reverse(config, unvec!, return_activity, tape, destination, x)
    function _shadow(annotation)

Reconstruct controls and flatten their cotangents into optimizer order.

# Arguments
- `config`, `destination`, `x`: configuration, control, and optimizer vector
- `tape`, `annotation`: reverse tape and duplicated annotation
# Output
- `result`: augmented result, reverse placeholders, or annotation shadow
"""
function augmented_primal(
    ::RevConfigWidth{1},
    func::Const{typeof(unvec!)},
    ::Type{<:Const},
    destination::Annotation{<:NamedTuple},
    x::Annotation{<:Vector},
)
    func.val(destination.val, x.val)
    return AugmentedReturn(nothing, nothing, nothing)
end
_shadow(annotation::Duplicated) = annotation.dval
_shadow(annotation::MixedDuplicated) = annotation.dval[]

function reverse(
    ::RevConfigWidth{1}, ::Const{typeof(unvec!)}, ::Type{<:Const}, ::Any,
    destination::Annotation, x::Annotation,
)
    destination isa Const && return (nothing, nothing)
    x isa Const || (_shadow(x) .+= vec(_shadow(destination)))
    # `unvec!` overwrites the wet values, so their cotangent is consumed. Zero only
    # those entries: with runtime activity, inactive shadow fields can alias primal
    # grid arrays, so `remake_zero!` on the whole struct is unsafe.
    unvec!(_shadow(destination), zero(vec(_shadow(destination))))
    return (nothing, nothing)
end

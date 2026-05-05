"""
    augmented_primal(config, unvec, Duplicated, template, uvec)

Forward pass for `y = unvec(template::BoundaryCondition, uvec)`.

`unvec` is the inverse layout operation to `vec`: it rebuilds a
`BoundaryCondition` from the compact vector of wet boundary values. The forward
rule returns `y` and allocates `gy = dJ/dy` for the reverse pass.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    y   = needs_primal(config) ? func.val(template.val, uvec.val)            : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

"""
    reverse(config, unvec, Duplicated, gy, template, uvec)

Reverse pass for boundary-condition `unvec`.

The adjoint of `unvec` is `vec`: `gy = dJ/dy` is read from the wet boundary
entries and accumulated into `guvec = dJ/duvec`, stored as `uvec.dval`.
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

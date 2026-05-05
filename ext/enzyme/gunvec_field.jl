"""
    augmented_primal(config, unvec, Duplicated, template, uvec)

Forward pass for `y = unvec(template::Field, uvec)`.

`unvec` is the inverse to `vec`. It rebuilds a `Field` from 
uvec. The rule returns `y` and allocates `gy = dJ/dy` for the reverse pass.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{Field{T,R,N,F}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, F <: AbstractArray{T,N}}
    y   = needs_primal(config) ? func.val(template.val, uvec.val)            : nothing
    gy = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, gy, gy)
end

"""
    reverse(config, unvec, Duplicated, gy, template, uvec)

Reverse pass for field `unvec`.

The adjoint of `unvec` is `vec`: `gy = dJ/dy` is read from wet grid cells and
accumulated into `guvec = dJ/duvec`, stored as `uvec.dval`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    gy::Field{T,R,N,F},
    template::Const{Field{T,R,N,F}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, F <: AbstractArray{T,N}}
    uvec.dval .+= gy.tracer[template.val.γ.wet]
    return (nothing, nothing)
end

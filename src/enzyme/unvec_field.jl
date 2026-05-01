"""
    augmented_primal(::Const{typeof(unvec)}, ..., ::Const{Field}, uvec)

Custom forward rule for `y = unvec(template::Field, uvec)`.

This is a scatter operation. See the tutorial for a full explanation.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{Field{T,R,N,F}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, F <: AbstractArray{T,N}}
    y   = needs_primal(config) ? func.val(template.val, uvec.val)            : nothing
    g_y = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, g_y, g_y)
end

"""
    reverse(::Const{typeof(unvec)}, ..., g_y, ::Const{Field}, uvec)

Custom reverse rule for `y = unvec(template::Field, uvec)`.

This is a gather operation, the adjoint of the forward scatter. See the tutorial.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    g_y::Field{T,R,N,F},
    template::Const{Field{T,R,N,F}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, F <: AbstractArray{T,N}}
    uvec.dval .+= g_y.tracer[template.val.γ.wet]
    return (nothing, nothing)
end

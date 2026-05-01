"""
    augmented_primal(::Const{typeof(unvec)}, ..., ::Const{BoundaryCondition}, uvec)

Custom forward rule for `y = unvec(template::BoundaryCondition, uvec)`.

This is a scatter operation. See the tutorial for a full explanation.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    y   = needs_primal(config) ? func.val(template.val, uvec.val)            : nothing
    g_y = needs_shadow(config) ? Enzyme.make_zero(func.val(template.val, uvec.val)) : nothing
    return AugmentedReturn(y, g_y, g_y)
end

"""
    reverse(::Const{typeof(unvec)}, ..., g_y, ::Const{BoundaryCondition}, uvec)

Custom reverse rule for `y = unvec(template::BoundaryCondition, uvec)`.

This is a gather operation, the adjoint of the forward scatter. See the tutorial.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(unvec)},
    ::Type{<:Duplicated},
    g_y::BoundaryCondition{T,R,N,G,B},
    template::Const{BoundaryCondition{T,R,N,G,B}},
    uvec::Duplicated{Vector{T}},
) where {T <: Real, R <: Real, N, G <: Integer, B <: AbstractArray{T,N}}
    if ndims(g_y.tracer) == 0
        # Handle 0-d boundary case (e.g., single surface cell in a 1-D grid).
        uvec.dval[1] += g_y.tracer[]
    else
        uvec.dval .+= g_y.tracer[template.val.wet]
    end
    return (nothing, nothing)
end

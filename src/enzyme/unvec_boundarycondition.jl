"""
    function augmented_primal(config, func, ::Type{<:Duplicated}, template, uvec)::AugmentedReturn

      Custom forward rule for `y = unvec(template::BoundaryCondition, uvec)`.

      This method is the forward pass of a scatter operation. It creates the
      full `BoundaryCondition` object from a vector of active values. The
      mathematical pattern for this rule (scatter/gather) is explained in the tutorial.

# Arguments
- `config`::RevConfigWidth{1}: Enzyme's internal configuration. `needs_primal` and
  `needs_shadow` query this to manage allocations.
- `func`::Const{typeof(unvec)}: The `unvec` function, wrapped in `Enzyme.Const` to
  mark it as a non-differentiable constant.
- `::Type{<:Duplicated}`: An unnamed type argument from Enzyme specifying that the
  function's output is differentiable.
- `template`::Const{BoundaryCondition}: The template object, which is not differentiated.
- `uvec`::Duplicated{Vector}: The vector of active values to be scattered, marked as
  differentiable. Its gradient will be stored in `uvec.dval`.

# Output
- `y`::AugmentedReturn: An `Enzyme.AugmentedReturn` struct containing the primal
  output `y` and its gradient storage `g_y`.
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
    function reverse(config, func, ::Type{<:Duplicated}, g_y, template, uvec)::(Nothing, Nothing)

      Custom reverse rule for `y = unvec(template::BoundaryCondition, uvec)`.

      This method performs a gather operation, which is the adjoint of the forward
      scatter. It accumulates gradients from the full field `g_y` back into the
      vector of active values `uvec.dval`.

      The implementation includes a special case for 1-D grids where the boundary
      condition tracer can be a 0-dimensional array (a single scalar). The `if` statement
      handles this by directly accessing the scalar gradient via `g_y.tracer[]`, while
      the `else` block handles the general N-dimensional case by using the `wet` mask
      to gather gradients from all active boundary points.

# Arguments
- `config`::RevConfigWidth{1}: Enzyme's internal configuration object.
- `func`::Const{typeof(unvec)}: The original function, marked as a constant.
- `::Type{<:Duplicated}`: An unnamed type argument from Enzyme specifying the
  differentiability of the original function's output.
- `g_y`::BoundaryCondition: Adjoint of the output from the forward pass. By the time
  this method is called, Enzyme has already accumulated the total gradient into `g_y`.
- `template`::Const{BoundaryCondition}: The non-differentiable template object. Its
  `wet` mask is used for indexing in the N-dimensional case.
- `uvec`::Duplicated{Vector}: The original vector of active values. Its `dval`
  field is updated in-place with the gathered gradients.

# Output
- `(nothing, nothing)`: Gradients are accumulated in-place.
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
        uvec.dval .+= g_y.tracer[template.val.wet] #this puts the gradietn in the correct vector form
    end
    return (nothing, nothing)
end

"""
    function augmented_primal(config, func, ::Type{<:Duplicated}, m, γ)::AugmentedReturn

      Custom forward rule for `A = watermassmatrix(m, γ)`.

      This method executes the forward pass, creating the sparse matrix `A` and
      allocating storage for its gradient. The mathematical pattern for this
      rule is explained in the tutorial.

# Arguments
- `config`::RevConfigWidth{1}: Enzyme's internal configuration object.
- `func`::Const{typeof(watermassmatrix)}: The function being differentiated, marked
  as a non-differentiable constant.
- `::Type{<:Duplicated}`: An unnamed type argument from Enzyme specifying that the
  function's output is differentiable.
- `m`::Duplicated: The collection of mass fractions, marked as differentiable.
- `γ`::Const{<:Grid}: The grid, which is not differentiated.

# Output
- `A`::AugmentedReturn: An `Enzyme.AugmentedReturn` struct containing the primal
  sparse matrix `A` and its gradient storage `g_A`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    m::Duplicated,
    γ::Const{<:Grid},
)
    A   = needs_primal(config) ? func.val(m.val, γ.val)            : nothing
    g_A = needs_shadow(config) ? Enzyme.make_zero(func.val(m.val, γ.val)) : nothing
    return AugmentedReturn(A, g_A, g_A)
end

"""
    function reverse(config, func, ::Type{<:Duplicated}, g_A, m, γ)::(Nothing, Nothing)

      Custom reverse rule for `A = watermassmatrix(m, γ)`.

      This method propagates gradients from the sparse matrix `g_A` backward,
      accumulating them into the `dval` fields of the mass fractions `m`.
      It does this by replaying the same loop structure as the forward pass
      to ensure the gradients from `g_A` are mapped to the correct mass fraction entries.

# Arguments
- `config`::RevConfigWidth{1}: Enzyme's internal configuration object.
- `func`::Const{typeof(watermassmatrix)}: The original function, marked as a constant.
- `::Type{<:Duplicated}`: An unnamed type argument from Enzyme specifying the
  differentiability of the original function's output.
- `g_A`: The adjoint of the sparse matrix from the forward pass, filled by
  subsequent rules in Enzyme's reverse pass.
- `m`::Duplicated: The mass fractions collection. Its `dval` field is updated
  in-place with the computed gradients.
- `γ`::Const{<:Grid}: The grid, which is not differentiated.

# Output
- `(nothing, nothing)`: Gradients are accumulated in-place.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    g_A,
    m::Duplicated,
    γ::Const{<:Grid},
)
    γval = γ.val
    R = γval.R
    Iint = cartesianindex(γval.interior)
    for I in Iint
        for (m1_val, m1_grad) in zip(m.val, m.dval)
            if m1_val.γ.wet[I]
                Istep, _ = step_cartesian(I, m1_val.position, γval)
                # Adjoint of: A[R[I], R[Istep]] = -m1.fraction[I]
                m1_grad.fraction[I] -= g_A[R[I], R[Istep]]
            end
        end
    end
    return (nothing, nothing)
end

"""
    augmented_primal(::Const{typeof(watermassmatrix)}, ..., m, γ)

Custom forward rule for `A = watermassmatrix(m, γ)`.

See the tutorial for a full explanation of the math.
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
    reverse(::Const{typeof(watermassmatrix)}, ..., g_A, m, γ)

Custom reverse rule for `A = watermassmatrix(m, γ)`.

This rule reads the gradient from `g_A` and accumulates it into the gradient
for the mass fractions in `m`. See the tutorial for a full derivation.
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

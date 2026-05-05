"""
    augmented_primal(config, watermassmatrix, Duplicated, m, γ)

Forward pass for `A = watermassmatrix(m, γ)`.

The mass fractions `m` are active. The grid `γ` is constant. The rule returns
the sparse matrix `A` and allocates `gA = dJ/dA` for the reverse pass.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    m::Duplicated,
    γ::Const{<:Grid},
)
    A   = needs_primal(config) ? func.val(m.val, γ.val)            : nothing
    gA = needs_shadow(config) ? Enzyme.make_zero(func.val(m.val, γ.val)) : nothing
    return AugmentedReturn(A, gA, gA)
end

"""
    reverse(config, watermassmatrix, Duplicated, gA, m, γ)

Reverse pass for `watermassmatrix`.

Here `gA = dJ/dA`. The rule replays the sparse-matrix assembly pattern and
accumulates each matching entry into `gm = dJ/dm`, stored in `m.dval`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    gA,
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
                m1_grad.fraction[I] += (-gA[R[I], R[Istep]])
            end
        end
    end
    return (nothing, nothing)
end

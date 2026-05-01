"""
    augmented_primal(::Const{typeof(\)}, ..., A, d)

Custom forward rule for the sparse linear solve `c = A \ d`.

See the tutorial for a full explanation of the math. This rule computes the
primal solution `c` and allocates storage for its gradient `g_c`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(\)},
    ::Type{<:Duplicated},
    A::Annotation{<:SparseMatrixCSC},
    d::Annotation{<:Field},
)
    c   = needs_primal(config) ? func.val(A.val, d.val)            : nothing
    g_c = needs_shadow(config) ? Enzyme.make_zero(func.val(A.val, d.val)) : nothing

    # Save c for the reverse pass and g_c for gradient accumulation.
    return AugmentedReturn(c, g_c, (c, g_c))
end

"""
    reverse(::Const{typeof(\)}, ..., tape, A, d)

Custom reverse rule for the sparse linear solve `c = A \ d`.

This implements the adjoint method. See the tutorial for a full derivation.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    tape,
    A::Annotation{<:SparseMatrixCSC},
    d::Annotation{<:Field},
)
    c, g_c = tape
    γ = c.γ
    g_c_vec = g_c.tracer[γ.wet]

    # Solve the adjoint system: A' * λ = g_c.
    # The adjoint variable λ is also the gradient g_d.
    g_d_vec = A.val' \ g_c_vec

    if d isa Duplicated
        d.dval.tracer[γ.wet] .+= g_d_vec
    end

    if A isa Duplicated
        c_vec = c.tracer[γ.wet]
        rows = SparseArrays.rowvals(A.dval)
        nzv  = SparseArrays.nonzeros(A.dval)
        for j in 1:size(A.dval, 2)
            for ki in SparseArrays.nzrange(A.dval, j)
                i = rows[ki]
                # g_A[i,j] -= λ[i] * c[j]
                nzv[ki] -= g_d_vec[i] * c_vec[j]
            end
        end
    end

    return (nothing, nothing)
end

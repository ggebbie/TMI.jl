"""
    augmented_primal(config, \\, return_activity, Alu, d)
    reverse(config, \\, return_activity, tape, Alu, d)

Solve `A c = d` for a `Field` `d` with TMI's `\\`, and propagate the cotangents
of `d` and of the water-mass matrix, `Ā = -d̄ cᵀ`, which goes into the
factorization's shadow on the sparsity pattern of `A`.

# Arguments
- `config`: reverse configuration
- `Alu`: duplicated UMFPACK factorization of `A`
- `d`: right-hand side `Field`
- `tape`: solution `c` and its cotangent

# Output
- `result`: augmented solve or reverse placeholders
"""
function augmented_primal(
    config::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    Alu::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    d::Annotation{<:Field},
)
    c = Alu.val \ d.val
    gc = Enzyme.make_zero(c)
    return AugmentedReturn(needs_primal(config) ? c : nothing, needs_shadow(config) ? gc : nothing,
        (c, gc))
end

function reverse(
    ::RevConfigWidth{1}, ::Const{typeof(\)}, ::Type{<:Duplicated}, tape,
    Alu::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    d::Union{Const{<:Field},Duplicated{<:Field}},
)
    c, gc = tape
    wet = d.val.γ.wet
    gd = Alu.val' \ gc.tracer[wet]
    d isa Duplicated && (d.dval.tracer[wet] .+= gd)
    # Ā = -d̄ cᵀ on the sparsity pattern of A; UMFPACK stores 0-based indices
    F, gF, cwet = Alu.val, Alu.dval, c.tracer[wet]
    for column in 1:F.n, k in (F.colptr[column] + 1):F.colptr[column + 1]
        gF.nzval[k] -= gd[F.rowval[k] + 1] * cwet[column]
    end
    return (nothing, nothing)
end

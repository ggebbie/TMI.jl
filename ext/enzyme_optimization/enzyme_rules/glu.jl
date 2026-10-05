"""
    Enzyme.make_zero(Alu)
    augmented_primal(config, lu!, return_activity, F, A)
    reverse(config, lu!, return_activity, _, F, A)

Create and propagate the sparse cotangent shadow for the UMFPACK refactorization
`lu!(F, A)`, which reuses the symbolic analysis (and ordering) of `F`. When `F`
already holds the values of `A`, as when Ipopt asks for the gradient at the
point whose cost it just evaluated, the refactorization is skipped.

# Arguments
- `Alu`, `F`: UMFPACK factorization and its duplicated shadow
- `A`: duplicated sparse source matrix
- `config`: reverse configuration

# Output
- `result`: factorization shadow, `AugmentedReturn`, or reverse placeholder
"""
function Enzyme.make_zero(Alu::SparseArrays.UMFPACK.UmfpackLU{Tv, Ti}) where {Tv, Ti}
    return typeof(Alu)(
        Alu.symbolic,
        Alu.numeric,
        Alu.m,
        Alu.n,
        copy(Alu.colptr),
        copy(Alu.rowval),
        zeros(Tv, length(Alu.nzval)),
        Alu.status,
        Alu.workspace,
        copy(Alu.control),
        copy(Alu.info),
        ReentrantLock(),
    )
end
function augmented_primal(
    config::RevConfigWidth{1},
    ::Const{typeof(lu!)},
    ::Type{<:Union{Const,Duplicated,Enzyme.DuplicatedNoNeed}},
    F::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    A::Duplicated{<:SparseMatrixCSC},
)
    nonzeros(A.val) == F.val.nzval || lu!(F.val, A.val)
    primal = needs_primal(config) ? F.val : nothing
    shadow = needs_shadow(config) ? F.dval : nothing
    return AugmentedReturn(primal, shadow, nothing)
end
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(lu!)},
    ::Type{<:Union{Const,Duplicated,Enzyme.DuplicatedNoNeed}},
    _,
    F::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    A::Duplicated{<:SparseMatrixCSC},
)
    # lu! overwrites F, so its cotangent passes entirely to A
    A.dval.nzval .+= F.dval.nzval
    fill!(F.dval.nzval, 0)
    return (nothing, nothing)
end

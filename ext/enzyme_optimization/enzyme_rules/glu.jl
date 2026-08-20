"""
    Enzyme.make_zero(Alu)
    augmented_primal(config, lu, return_activity, A)
    reverse(config, lu, return_activity, gAlu, A)

Create and propagate the sparse cotangent shadow for UMFPACK factorization.

# Arguments
- `Alu`, `A`: factorization and duplicated sparse source matrix
- `config`, `gAlu`: reverse configuration and factorization cotangent

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
    func::Const{typeof(lu)},
    ::Type{<:Union{Duplicated,Enzyme.DuplicatedNoNeed}},
    A::Duplicated{<:SparseMatrixCSC},
)
    Alu = func.val(A.val)
    gAlu = needs_shadow(config) ? Enzyme.make_zero(Alu) : nothing
    primal = needs_primal(config) ? Alu : nothing
    return AugmentedReturn(primal, gAlu, gAlu)
end
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(lu)},
    ::Type{<:Union{Duplicated,Enzyme.DuplicatedNoNeed}},
    gAlu::SparseArrays.UMFPACK.UmfpackLU,
    A::Duplicated{<:SparseMatrixCSC},
)
    A.dval.nzval .+= gAlu.nzval
    return (nothing,)
end
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(lu)},
    ::Type{<:Union{Duplicated,Enzyme.DuplicatedNoNeed}},
    ::Nothing,
    ::Duplicated{<:SparseMatrixCSC},
)
    return (nothing,)
end

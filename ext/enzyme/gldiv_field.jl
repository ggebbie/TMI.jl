"""
    augmented_primal(config, left division, Duplicated, A, d)

Forward pass for solving a linear system with a `Field` right-hand side and
constant `A`.

Enzyme's reverse mode runs in two phases. The forward rule returns the primal
output `c`, allocates `gc = dJ/dc` when a reverse pass is needed, and tapes the
values needed by `reverse`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(\)},
    ::Type{<:Duplicated},
    A::Const{<:Union{SparseMatrixCSC, SparseArrays.UMFPACK.UmfpackLU}},
    d::Annotation{<:Field},
)
    γ = d.val.γ
    dvec = d.val.tracer[γ.wet]
    # Keep this on the wet-vector solve path. Calling `func.val(A.val, d.val)`
    # dispatches to `\(A, d::Field)`, which has been unstable in long Enzyme
    # optimization loops.
    cvec = A.val \ dvec
    c = zero(d.val)
    c.tracer[γ.wet] .= cvec
    primal = needs_primal(config) ? c : nothing
    gc = needs_shadow(config) ? Enzyme.make_zero(c) : nothing
    return AugmentedReturn(primal, gc, gc)
end

"""
    augmented_primal(config, left division, Duplicated, A, d)

Forward pass for solving a linear system with active `A`.

This method accepts either `A::SparseMatrixCSC` or `A::UmfpackLU`.  In the LU
case, `A.dval.nzval` is used as the sparse adjoint buffer on the original
matrix pattern stored by the factorization.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(\)},
    ::Type{<:Duplicated},
    A::Duplicated{<:Union{SparseMatrixCSC, SparseArrays.UMFPACK.UmfpackLU}},
    d::Annotation{<:Field},
)
    γ = d.val.γ
    dvec = d.val.tracer[γ.wet]
    # Same reason as above: avoid `\(A, d::Field)` dispatch in Enzyme loops.
    cvec = A.val \ dvec
    c = zero(d.val)
    c.tracer[γ.wet] .= cvec
    primal = needs_primal(config) ? c : nothing
    gc = needs_shadow(config) ? Enzyme.make_zero(c) : nothing
    return AugmentedReturn(primal, gc, (c, gc))
end

"""
    reverse(config, left division, Duplicated, gc, A, d)

Reverse pass for a linear solve with constant `A`.

Given `gc = dJ/dc`, the adjoint equation is `A' * gd = gc`, so
`gd` is computed by solving the transpose system.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    gc,
    A::Const{<:Union{SparseMatrixCSC, SparseArrays.UMFPACK.UmfpackLU}},
    d::Annotation{<:Field},
)
    γ = d.val.γ
    gcvec = gc.tracer[γ.wet]

    gdvec = A.val' \ gcvec

    if d isa Duplicated
        d.dval.tracer[γ.wet] .+= gdvec
    end

    return (nothing, nothing)
end

"""
    reverse(config, left division, Duplicated, tape, A, d)

Reverse pass for a linear solve with active `A`.

Accumulates `dJ/dd = A' \\ gc` and `dJ/dA = -gd*c'`.  For a sparse matrix the
adjoint is written to `nonzeros(A.dval)`; for an LU factorization it is written
to `A.dval.nzval`, which has the original matrix sparsity pattern.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    tape,
    A::Duplicated{<:Union{SparseMatrixCSC, SparseArrays.UMFPACK.UmfpackLU}},
    d::Annotation{<:Field},
)
    c, gc = tape
    γ = d.val.γ
    gcvec = gc.tracer[γ.wet]
    gdvec = A.val' \ gcvec

    if d isa Duplicated
        d.dval.tracer[γ.wet] .+= gdvec
    end

    cvec = c.tracer[γ.wet]
    if A.val isa SparseMatrixCSC
        Asparse = A.val
        nzv = SparseArrays.nonzeros(A.dval)
    elseif A.val isa SparseArrays.UMFPACK.UmfpackLU
        Asparse = SparseMatrixCSC(
            A.val.m,
            A.val.n,
            copy(A.val.colptr .+ 1),
            copy(A.val.rowval .+ 1),
            zeros(eltype(A.val.nzval), length(A.val.nzval)),
        )
        nzv = A.dval.nzval
    else
        error("unsupported active solve type: ", typeof(A.val))
    end

    rows = SparseArrays.rowvals(Asparse)
    for j in 1:size(Asparse, 2)
        for ki in SparseArrays.nzrange(Asparse, j)
            i = rows[ki]
            nzv[ki] -= gdvec[i] * cvec[j]
        end
    end

    return (nothing, nothing)
end

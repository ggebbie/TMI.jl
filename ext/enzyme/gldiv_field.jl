"""
    Enzyme.make_zero(A::SparseArrays.UMFPACK.UmfpackLU)

Create a zero shadow for an LU factorization. The differentiable payload for
the original matrix pattern is `nzval`; symbolic/numeric factorization handles
are metadata for the solve and are reused.
"""
function Enzyme.make_zero(A::SparseArrays.UMFPACK.UmfpackLU{Tv, Ti}) where {Tv, Ti}
    return typeof(A)(
        A.symbolic,
        A.numeric,
        A.m,
        A.n,
        copy(A.colptr),
        copy(A.rowval),
        zeros(Tv, length(A.nzval)),
        A.status,
        A.workspace,
        copy(A.control),
        copy(A.info),
        ReentrantLock(),
    )
end

"""
    `augmented_primal(config, \\, Duplicated, A, d)`

Custom Enzyme forward rule for TMI's left-division method:

    \\(A, d::Field)::Field

This method computes `c = A \\ d` with `d::Field` and constant `A`.

Enzyme's reverse mode runs in two phases. The forward rule returns the primal
output `c`, allocates `gc = dJ/dc` when a reverse pass is needed, and tapes the
values needed by `reverse`.

Linked TMI function:
- `TMI.\\(A, d::Field)`

Arguments:
- `A`: constant sparse/LU operator.
- `d`: active or mixed-active `Field` right-hand side.
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
    `augmented_primal(config, \\, Duplicated, A, d)`

Forward rule for the same method

    \\(A, d::Field)::Field

when `A` is active.

This method is for `A::Duplicated{SparseMatrixCSC}` and tapes `(c, gc)` so the
reverse pass can accumulate both `dJ/dd` and `dJ/dA`.

Linked TMI function:
- `TMI.\\(A, d::Field)`

Arguments:
- `A`: active sparse/LU operator.
- `d`: active or mixed-active `Field` right-hand side.
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
    `reverse(config, \\, Duplicated, gc, A, d)`

Reverse pass for

    \\(A, d::Field)::Field

with constant `A`.

Given `gc = dJ/dc`, the adjoint equation is `A' * gd = gc`, so
`gd = dJ/dd = A' \\ gc`.

Linked TMI function:
- `TMI.\\(A, d::Field)`

Arguments:
- `gc`: gradient w.r.t. solve output `c`.
- `A`: constant sparse/LU operator.
- `d`: input `Field` receiving `dJ/dd`.
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
    `reverse(config, \\, Duplicated, tape, A, d)`

Reverse pass for

    \\(A, d::Field)::Field

with active sparse `A`.

Accumulates `dJ/dd = A' \\ gc` and `dJ/dA = -gd*c'` on the stored sparse
pattern of `A.dval`.

Linked TMI function:
- `TMI.\\(A, d::Field)`

Arguments:
- `tape`: `(c, gc)` from the forward rule.
- `A`: active `SparseMatrixCSC` operator receiving `dJ/dA`.
- `d`: input `Field` receiving `dJ/dd`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    tape,
    A::Duplicated{<:SparseMatrixCSC},
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
    rows = SparseArrays.rowvals(A.val)
    nzv = SparseArrays.nonzeros(A.dval)
    for j in 1:size(A.val, 2)
        for ki in SparseArrays.nzrange(A.val, j)
            i = rows[ki]
            nzv[ki] -= gdvec[i] * cvec[j]
        end
    end

    return (nothing, nothing)
end

"""
    `reverse(config, \\, Duplicated, tape, A, d)` for `A::Duplicated{UmfpackLU}`

Same adjoint as the sparse-matrix case, but writes `dJ/dA` into `A.dval.nzval`
using UMFPACK's stored sparsity pattern.

Linked TMI function:
- `TMI.\\(A, d::Field)`

Arguments:
- `tape`: `(c, gc)` from the forward rule.
- `A`: active `UmfpackLU` operator receiving `dJ/dA`.
- `d`: input `Field` receiving `dJ/dd`.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    tape,
    A::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
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
    nzv = A.dval.nzval
    Asparse = SparseMatrixCSC(
        A.val.m,
        A.val.n,
        deepcopy(A.val.colptr .+ 1),
        deepcopy(A.val.rowval .+ 1),
        zeros(eltype(A.val.nzval), length(A.val.nzval)),
    )
    rows = SparseArrays.rowvals(Asparse)
    for j in 1:size(Asparse, 2)
        for ki in SparseArrays.nzrange(Asparse, j)
            i = rows[ki]
            nzv[ki] -= gdvec[i] * cvec[j]
        end
    end

    return (nothing, nothing)
end

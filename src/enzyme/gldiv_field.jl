const SparseOrLU = Union{SparseMatrixCSC, SparseArrays.UMFPACK.UmfpackLU}

"""
    augmented_primal(config, \\ , Duplicated, A, d)

Forward pass for `c = A \\ d`, where `d` is a `Field` and `A` is either a
`SparseMatrixCSC` or a UMFPACK LU factorization.

Enzyme's reverse mode runs in two phases. The forward rule returns the primal
output `c`, allocates `gc = dJ/dc` when a reverse pass is needed, and tapes the
values needed by `reverse`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(\)},
    ::Type{<:Duplicated},
    A::Annotation{<:SparseOrLU},
    d::Annotation{<:Field},
)
    c = func.val(A.val, d.val)
    primal = needs_primal(config) ? c : nothing
    gc = needs_shadow(config) ? Enzyme.make_zero(c) : nothing
    Asparse = needs_shadow(config) ? if A.val isa SparseMatrixCSC
            A.val
        else
            SparseMatrixCSC(A.val.m, A.val.n, A.val.colptr .+ 1, A.val.rowval .+ 1, copy(A.val.nzval))
        end : nothing
    return AugmentedReturn(primal, gc, (c, gc, Asparse))
end

"""
    reverse(config, \\ , Duplicated, tape, A, d)

Reverse pass for `c = A \\ d`.

Given `gc = dJ/dc`, the adjoint equation is `A' * gd = gc`, so
`gd = dJ/dd = A' \\ gc`. If `A` itself is active and its shadow is a sparse
matrix, the same solve gives `gA = dJ/dA = -gd * c'`, accumulated only on the
stored sparse pattern.
"""
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    tape,
    A::Annotation{<:SparseOrLU},
    d::Annotation{<:Field},
)
    c, gc, Asparse = tape
    γ = c.γ
    gcvec = gc.tracer[γ.wet]

    gdvec = A.val' \ gcvec

    if d isa Duplicated
        d.dval.tracer[γ.wet] .+= gdvec
    end

    if A isa Duplicated && A.dval isa SparseMatrixCSC
        cvec = c.tracer[γ.wet]
        rows = SparseArrays.rowvals(Asparse)
        nzv  = SparseArrays.nonzeros(A.dval)
        for j in 1:size(Asparse, 2)
            for ki in SparseArrays.nzrange(Asparse, j)
                i = rows[ki]
                nzv[ki] -= gdvec[i] * cvec[j]
            end
        end
    end

    return (nothing, nothing)
end

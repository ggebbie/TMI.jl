#=
Enzyme rules: the adjoints (TMI's `g` functions) of the four steps of the forward
model that Enzyme is not left to differentiate itself. Nobody calls these
directly; they matter only when one of these four functions changes. Enzyme
differentiates everything else in the cost function, including `observe`.
Each rule is a pair of methods that Enzyme calls in place of differentiating the
function: `augmented_primal` runs the function and keeps what the backward pass
needs, and `reverse` passes the cotangents back to the inputs.
In the order the forward model calls them:
- `unvec!` writes `uvec` into the controls `u`. Its rule passes their cotangents
  back to `uvec`, and zeros only the wet values, because with runtime activity
  the shadows of fixed controls can share memory with the first guesses.
- `watermassmatrix` builds the sparse water-mass matrix `A` from the mass
  fractions. Its rule takes the cotangent of `A` back to each fraction.
- `lu!` refactors `A` with UMFPACK, a C library Enzyme cannot differentiate. Its
  rule passes the factorization's cotangent on to `A`, and skips the
  refactorization when `F` already holds `A`.
- `\` solves `A c = d` with the UMFPACK factors. Its rule solves the adjoint
  system `Aᵀ λ = c̄` with the same factors.
=#

"""
    function EnzymeRules.augmented_primal(config, unvec!, return_activity, destination, x)
    function EnzymeRules.reverse(config, unvec!, return_activity, tape, destination, x)
    function _shadow(annotation)

Reconstruct controls and flatten their cotangents into optimizer order:
`ūvec += vec(ū)`, then `ū = 0` at the wet points.

# Arguments
- `config`, `destination`, `x`: configuration, control, and optimizer vector
- `tape`, `annotation`: reverse tape and duplicated annotation
# Output
- `result`: augmented result, reverse placeholders, or annotation shadow
"""
function EnzymeRules.augmented_primal(
    ::RevConfigWidth{1},
    func::Const{typeof(unvec!)},
    ::Type{<:Const},
    destination::Annotation{<:NamedTuple},
    x::Annotation{<:Vector},
)
    func.val(destination.val, x.val)
    return AugmentedReturn(nothing, nothing, nothing)
end
_shadow(annotation::Duplicated) = annotation.dval
_shadow(annotation::MixedDuplicated) = annotation.dval[]

function EnzymeRules.reverse(
    ::RevConfigWidth{1}, ::Const{typeof(unvec!)}, ::Type{<:Const}, ::Any,
    destination::Annotation, x::Annotation,
)
    destination isa Const && return (nothing, nothing)
    if !(x isa Const)
        _shadow(x) .+= vec(_shadow(destination))
    end
    # `unvec!` overwrites the wet values, so their cotangent is consumed. Zero only
    # those entries: with runtime activity, inactive shadow fields can alias primal
    # grid arrays, so `remake_zero!` on the whole struct is unsafe.
    unvec!(_shadow(destination), zero(vec(_shadow(destination))))
    return (nothing, nothing)
end

"""
    function EnzymeRules.augmented_primal(config, watermassmatrix, return_activity, m, γ)
    function EnzymeRules.reverse(config, watermassmatrix, return_activity, gA, m, γ)

Build the active water-mass matrix and reverse its cotangent to fractions:
`m̄ = -Ā` at each fraction's entry of `A`, which holds `-m` off the diagonal.

# Arguments
- `config`, `m`, `γ`: reverse configuration, duplicated fractions, and grid
- `gA`: sparse water-mass-matrix cotangent
# Output
- `result`: augmented matrix result or reverse placeholders
"""
function EnzymeRules.augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    m::Union{Const,Duplicated},
    γ::Annotation{<:Grid},
)
    A = func.val(m.val, γ.val)
    gA = Enzyme.make_zero(A)
    return AugmentedReturn(needs_primal(config) ? A : nothing, needs_shadow(config) ? gA : nothing, gA)
end
function EnzymeRules.reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    gA,
    m::Union{Const,Duplicated},
    γ::Annotation{<:Grid},
)
    # fixed mass fractions: A is constant
    m isa Const && return (nothing, nothing)
    grid = γ.val
    R = grid.R
    # fractions outer: taking an element of `m` in the inner loop allocates on Julia 1.13
    for (mₖ, gmₖ) in zip(m.val, m.dval)
        # with runtime activity, a fixed fraction arrives as its own shadow: zero derivative
        gmₖ.fraction === mₖ.fraction && continue
        for cell in cartesianindex(wet(mₖ))
            neighbor, _ = step_cartesian(cell, mₖ.position, grid)
            gmₖ.fraction[cell] -= gA[R[cell], R[neighbor]]
        end
    end
    return (nothing, nothing)
end

"""
    function Enzyme.make_zero(Alu)
    function EnzymeRules.augmented_primal(config, lu!, return_activity, F, A)
    function EnzymeRules.reverse(config, lu!, return_activity, _, F, A)

Create and propagate the sparse cotangent shadow for the UMFPACK refactorization
`lu!(F, A)`, which reuses the symbolic analysis (and ordering) of `F`. When `F`
already holds the values of `A`, as when Ipopt asks for the gradient at the
point whose cost it just evaluated, the refactorization is skipped. The reverse
pass is `Ā += F̄`, since `lu!` overwrites `F`.

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
function EnzymeRules.augmented_primal(
    config::RevConfigWidth{1},
    ::Const{typeof(lu!)},
    ::Type{<:Union{Const,Duplicated,Enzyme.DuplicatedNoNeed}},
    F::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    A::Duplicated{<:SparseMatrixCSC},
)
    # skip the refactorization when F already holds A
    if nonzeros(A.val) != F.val.nzval
        lu!(F.val, A.val)
    end
    primal = needs_primal(config) ? F.val : nothing
    shadow = needs_shadow(config) ? F.dval : nothing
    return AugmentedReturn(primal, shadow, nothing)
end
function EnzymeRules.reverse(
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

"""
    function EnzymeRules.augmented_primal(config, \\, return_activity, Alu, d)
    function EnzymeRules.reverse(config, \\, return_activity, tape, Alu, d)

Solve `A c = d` for a `Field` `d` with TMI's `\\`, and propagate the cotangents:
`λ = A⁻ᵀ c̄`, `d̄ += λ`, and `Ā -= λ cᵀ` on the sparsity pattern of `A`, which goes
into the factorization's shadow.

# Arguments
- `config`: reverse configuration
- `Alu`: duplicated UMFPACK factorization of `A`
- `d`: right-hand side `Field`
- `tape`: solution `c` and its cotangent
# Output
- `result`: augmented solve or reverse placeholders
"""
function EnzymeRules.augmented_primal(
    config::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    Alu::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    d::Union{Const{<:Field},Duplicated{<:Field}},
)
    c = Alu.val \ d.val
    gc = Enzyme.make_zero(c)
    return AugmentedReturn(needs_primal(config) ? c : nothing, needs_shadow(config) ? gc : nothing,
        (c, gc))
end

function EnzymeRules.reverse(
    ::RevConfigWidth{1}, ::Const{typeof(\)}, ::Type{<:Duplicated}, tape,
    Alu::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    d::Union{Const{<:Field},Duplicated{<:Field}},
)
    c, gc = tape
    wet = d.val.γ.wet
    gd = Alu.val' \ gc.tracer[wet]
    if d isa Duplicated
        d.dval.tracer[wet] .+= gd
    end
    # Ā = -d̄ cᵀ on the sparsity pattern of A; UMFPACK stores 0-based indices
    F, gF, cwet = Alu.val, Alu.dval, c.tracer[wet]
    for column in 1:F.n, k in (F.colptr[column] + 1):F.colptr[column + 1]
        gF.nzval[k] -= gd[F.rowval[k] + 1] * cwet[column]
    end
    return (nothing, nothing)
end

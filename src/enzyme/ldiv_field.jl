"""
    function augmented_primal(config, func, ::Type{<:Duplicated}, A, d)::AugmentedReturn

      Custom forward rule for the sparse linear solve `c = A \\ d`.

      This method is called by Enzyme's reverse mode to execute the forward pass.
      It computes the primal solution `c` and prepares the necessary storage for the
      reverse pass. The mathematical pattern for this rule (the adjoint method)
      is explained in the tutorial.

# Arguments
- `config`::RevConfigWidth{1}: Enzyme's internal configuration for the rule. The `needs_primal`
  and `needs_shadow` functions query this object to determine if the primal value or
  gradient storage are required for the current differentiation task.
- `func`::Const{typeof(\\)}: The function being differentiated (``), wrapped in `Enzyme.Const`
  to mark it as a non-differentiable constant.
- `::Type{<:Duplicated}`: An unnamed type argument from Enzyme specifying that the
  function's output is differentiable.
- `A`::Annotation{<:SparseMatrixCSC}: The sparse matrix argument, wrapped in
  `Enzyme.Annotation` which provides metadata about the argument's differentiability.
- `d`::Annotation{<:Field}: The right-hand-side field argument, also wrapped in
  `Enzyme.Annotation`.

# Output
- `c`::AugmentedReturn: An `Enzyme.AugmentedReturn` struct that contains the primal
  output (`c`), storage for the output's gradient (`g_c`), and any values that need
  to be "taped" for the reverse pass.
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
    function reverse(config, func, ::Type{<:Duplicated}, tape, A, d)::(Nothing, Nothing)

      Custom reverse rule for the sparse linear solve `c = A \\ d`.

      This method implements the adjoint method to propagate gradients backward.
      It first solves the adjoint system `A' * λ = g_c` to get the adjoint variable
      `λ` (which is also the gradient `g_d`). It then uses `λ` to accumulate the
      gradients for `A` and `d`. The `if` blocks ensure that gradients are only
      computed and stored for arguments that were marked as active (`Duplicated`).

# Arguments
- `config`::RevConfigWidth{1}: Enzyme's internal configuration object.
- `func`::Const{typeof(\\)}: The original function, marked as a constant.
- `::Type{<:Duplicated}`: An unnamed type argument from Enzyme specifying the
  differentiability of the original function's output.
- `tape`: Contains the values saved from `AugmentedReturn` in the forward pass,
  in this case the primal solution `c` and its gradient storage `g_c`.
- `A`::Annotation{<:SparseMatrixCSC}: The sparse matrix argument. If it is a
  `Duplicated` annotation, its `dval` field will be updated with the gradient.
- `d`::Annotation{<:Field}: The right-hand-side field argument. If it is a
  `Duplicated` annotation, its `dval` field will be updated with the gradient.

# Output
- `(nothing, nothing)`: This rule returns no gradients directly. Gradients are
  accumulated in-place into the `dval` fields of the `Duplicated` arguments.
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

    # Solve the adjoint system: A' * g_d = g_c.
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

"""
    Custom reverse-mode rule for `Alu \\ d::Field` with a constant UMFPACK LU factorization.

The factorization is treated as non-differentiable, so the reverse pass only propagates
the adjoint solve into the right-hand side: `g_d = Alu' \\ g_c`.
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(\)},
    ::Type{<:Duplicated},
    Alu::Annotation{<:SparseArrays.UMFPACK.UmfpackLU},
    d::Annotation{<:Field},
)
    dval = d.val
    γ = dval.γ
    tracer = copy(dval.tracer)
    tracer[γ.wet] .= Alu.val \ dval.tracer[γ.wet]
    c = Field(tracer, γ, dval.name, dval.longname, dval.units)
    primal = needs_primal(config) ? c : nothing
    g_c = needs_shadow(config) ? Enzyme.make_zero(c) : nothing

    return AugmentedReturn(primal, g_c, g_c)
end

function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    g_c,
    Alu::Annotation{<:SparseArrays.UMFPACK.UmfpackLU},
    d::Annotation{<:Field},
)
    if d isa Duplicated
        γ = d.val.γ
        g_c_vec = g_c.tracer[γ.wet]
        g_d_vec = Alu.val' \ g_c_vec
        d.dval.tracer[γ.wet] .+= g_d_vec
    end

    return (nothing, nothing)
end

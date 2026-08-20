"""
    augmented_primal(config, \\, return_activity, Alu, right_hand_side)
    reverse(config, \\, return_activity, tape, Alu, right_hand_side)
    _solvetape(Alu, c, gc)

Solve `A c = d` and propagate field and sparse-factorization cotangents.

# Arguments
- `config`, `Alu`, `right_hand_side`: configuration, factorization, and field
- `tape`, `c`, `gc`: reverse tape, solution, and cotangent

# Output
- `result`: augmented solve, typed tape, or reverse placeholder
"""
_solvetape(::Const, ::Field, gc) = gc
_solvetape(::Duplicated, c::Field, gc) = (c, gc)
function augmented_primal(
    config::RevConfigWidth{1},
    ::Const{typeof(\)},
    ::Type{<:Duplicated},
    Alu::Union{Const{<:SparseArrays.UMFPACK.UmfpackLU},
        Duplicated{<:SparseArrays.UMFPACK.UmfpackLU}},
    right_hand_side::Annotation{<:Field},
)
    wet_mask = right_hand_side.val.γ.wet
    c = zero(right_hand_side.val)
    c.tracer[wet_mask] .= Alu.val \ right_hand_side.val.tracer[wet_mask]
    primal = needs_primal(config) ? c : nothing
    gc = needs_shadow(config) ? Enzyme.make_zero(c) : nothing
    return AugmentedReturn(primal, gc, _solvetape(Alu, c, gc))
end
function reverse(
    ::RevConfigWidth{1}, ::Const{typeof(\)}, ::Type{<:Duplicated}, tape,
    Alu::Const{<:SparseArrays.UMFPACK.UmfpackLU},
    right_hand_side::Duplicated{<:Field},
)
    wet_mask = right_hand_side.val.γ.wet
    right_hand_side.dval.tracer[wet_mask] .+= Alu.val' \ tape.tracer[wet_mask]
    return (nothing, nothing)
end
function reverse(
    ::RevConfigWidth{1}, ::Const{typeof(\)}, ::Type{<:Duplicated}, ::Any,
    ::Const{<:SparseArrays.UMFPACK.UmfpackLU}, ::Const{<:Field},
)
    return (nothing, nothing)
end
function reverse(
    ::RevConfigWidth{1}, ::Const{typeof(\)}, ::Type{<:Duplicated}, tape,
    Alu::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    right_hand_side::Duplicated{<:Field},
)
    c, gc = tape
    wet_mask = right_hand_side.val.γ.wet
    gd = Alu.val' \ gc.tracer[wet_mask]
    right_hand_side.dval.tracer[wet_mask] .+= gd
    c_wet = c.tracer[wet_mask]
    for column in 1:Alu.val.n
        nonzero_range = (Alu.val.colptr[column] + 1):Alu.val.colptr[column + 1]
        for nonzero_index in nonzero_range
            row = Alu.val.rowval[nonzero_index] + 1
            Alu.dval.nzval[nonzero_index] -= gd[row] * c_wet[column]
        end
    end
    return (nothing, nothing)
end
function reverse(
    ::RevConfigWidth{1}, ::Const{typeof(\)}, ::Type{<:Duplicated}, tape,
    Alu::Duplicated{<:SparseArrays.UMFPACK.UmfpackLU},
    right_hand_side::Const{<:Field},
)
    c, gc = tape
    wet_mask = right_hand_side.val.γ.wet
    gd = Alu.val' \ gc.tracer[wet_mask]
    c_wet = c.tracer[wet_mask]
    for column in 1:Alu.val.n
        nonzero_range = (Alu.val.colptr[column] + 1):Alu.val.colptr[column + 1]
        for nonzero_index in nonzero_range
            row = Alu.val.rowval[nonzero_index] + 1
            Alu.dval.nzval[nonzero_index] -= gd[row] * c_wet[column]
        end
    end
    return (nothing, nothing)
end

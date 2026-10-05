"""
    augmented_primal(config, watermassmatrix, return_activity, m, γ)
    reverse(config, watermassmatrix, return_activity, gA, m, γ)

Build the active water-mass matrix and reverse its cotangent to fractions.

# Arguments
- `config`, `m`, `γ`: reverse configuration, duplicated fractions, and grid
- `gA`: sparse water-mass-matrix cotangent

# Output
- `result`: augmented matrix result or reverse placeholders
"""
function augmented_primal(
    config::RevConfigWidth{1},
    func::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    m::Union{Const,Duplicated},
    γ::Annotation{<:Grid},
)
    needs_matrix = needs_primal(config) || needs_shadow(config)
    A = needs_matrix ? func.val(m.val, γ.val) : nothing
    primal = needs_primal(config) ? A : nothing
    gA = needs_shadow(config) ? Enzyme.make_zero(A) : nothing
    return AugmentedReturn(primal, gA, gA)
end
function reverse(
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
        for cell in cartesianindex(wet(mₖ))
            neighbor, _ = step_cartesian(cell, mₖ.position, grid)
            gmₖ.fraction[cell] -= gA[R[cell], R[neighbor]]
        end
    end
    return (nothing, nothing)
end

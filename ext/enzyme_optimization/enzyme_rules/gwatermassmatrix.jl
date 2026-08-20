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
    mass_fractions::Duplicated,
    γ::Annotation{<:Grid},
)
    needs_matrix = needs_primal(config) || needs_shadow(config)
    A = needs_matrix ? func.val(mass_fractions.val, γ.val) : nothing
    primal = needs_primal(config) ? A : nothing
    gA = needs_shadow(config) ? Enzyme.make_zero(A) : nothing
    return AugmentedReturn(primal, gA, gA)
end
function reverse(
    ::RevConfigWidth{1},
    ::Const{typeof(watermassmatrix)},
    ::Type{<:Duplicated},
    gA,
    mass_fractions::Duplicated,
    γ::Annotation{<:Grid},
)
    grid = γ.val
    R = grid.R
    for cell in cartesianindex(grid.interior)
        for (fraction, g_fraction) in zip(mass_fractions.val, mass_fractions.dval)
            wet(fraction)[cell] || continue
            neighbor, _ = step_cartesian(cell, fraction.position, grid)
            g_fraction.fraction[cell] -= gA[R[cell], R[neighbor]]
        end
    end
    return (nothing, nothing)
end

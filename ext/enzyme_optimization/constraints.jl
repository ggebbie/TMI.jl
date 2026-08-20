"""
    mass_conservation_constraints!(constraint_values, m, template)
    mass_conservation_constraints(inversion)
    mass_conservation_constraints(inversion, control_offset, control_scale)

Evaluate and configure directional mass-conservation constraints.

# Arguments
- `constraint_values`, `m`, `template`: destination, values, and templates
- `inversion`: inversion definition
- `control_offset`, `control_scale`: optimizer recentering and normalization

# Output
- `result`: in-place result or callbacks, bounds, and Jacobian prototype
"""
function mass_conservation_constraints!(constraint_values, m, template)
    mass_fractions = unvec(template, m)
    interior_cells = cartesianindex(first(mass_fractions).γ.interior)
    for (row, cell) in enumerate(interior_cells)
        total_fraction = zero(eltype(m))
        for fraction in mass_fractions
            wet(fraction)[cell] || continue
            total_fraction += fraction.fraction[cell]
        end
        constraint_values[row] = total_fraction
    end
    return nothing
end

mass_conservation_constraints(inversion::Inversion) =
    mass_conservation_constraints(inversion, zero(inversion.x₀),
        ones(eltype(inversion.x₀), length(inversion.x₀)))

function mass_conservation_constraints(inversion::Inversion,
    control_offset::AbstractVector, control_scale::AbstractVector)
    if isnothing(inversion.control_priors.mass_fraction)
        return (callbacks=(;), bounds=(;))
    end
    m_range = inversion.control_ranges.mass_fraction
    mass_fraction_offset = control_offset[m_range]
    mass_fraction_scale = control_scale[m_range]
    template = inversion.control_templates.mass_fractions
    γ = first(template).γ
    row_at = zeros(Int, size(γ.interior))
    row_at[γ.interior] = eachindex(cartesianindex(γ.interior))
    jacobian_rows = Int[]
    for fraction in template
        for cell in cartesianindex(fraction.γ.wet)
            push!(jacobian_rows, row_at[cell])
        end
    end
    mass_fraction_jacobian = sparse(jacobian_rows, eachindex(jacobian_rows),
        mass_fraction_scale, sum(γ.interior), length(jacobian_rows))
    left_zeros = spzeros(eltype(mass_fraction_jacobian),
        size(mass_fraction_jacobian, 1), first(m_range) - 1)
    right_zeros = spzeros(eltype(mass_fraction_jacobian),
        size(mass_fraction_jacobian, 1), length(inversion.x₀) - last(m_range))
    jacobian_prototype = hcat(left_zeros, mass_fraction_jacobian, right_zeros)
    physical_mass_fraction_controls = similar(mass_fraction_offset)
    constraint! = (values, optimizer_controls, _) -> begin
        @views @. physical_mass_fraction_controls = mass_fraction_offset +
            mass_fraction_scale * optimizer_controls[m_range]
        mass_conservation_constraints!(values,
            physical_mass_fraction_controls, template)
    end
    jacobian! = (jacobian, _, _) -> begin
        copyto!(nonzeros(jacobian), nonzeros(jacobian_prototype))
        nothing
    end
    target = ones(eltype(inversion.x₀), size(jacobian_prototype, 1))
    callbacks = (cons=constraint!, cons_j=jacobian!,
        cons_jac_prototype=jacobian_prototype)
    return (; callbacks, bounds=(lcons=target, ucons=target))
end

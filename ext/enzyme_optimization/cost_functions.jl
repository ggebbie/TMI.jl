"""
    data_cost(r, Q_y)

Return the precision-weighted squared observation residual.

# Arguments
- `r`, `Q_y`: residual vector and observation precision matrix

# Output
- `J_y`: scalar data-misfit cost
"""
data_cost(r, Q_y) = dot(r, Q_y, r)

"""
    boundary_prior_cost(Δb, Q_b)

Return the precision-weighted squared boundary deviation.

# Arguments
- `Δb`, `Q_b`: boundary deviation and prior precision matrix

# Output
- `J_b`: scalar boundary-prior cost
"""
boundary_prior_cost(Δb, Q_b) = dot(Δb, Q_b, Δb)

"""
    boundary_smoothness_cost(Δb, Q_s)

Return the precision-weighted boundary roughness cost.

# Arguments
- `Δb`, `Q_s`: boundary deviation and smoothness precision matrix

# Output
- `J_s`: scalar boundary-smoothness cost
"""
boundary_smoothness_cost(Δb, Q_s) = dot(Δb, Q_s, Δb)

"""
    source_prior_cost(Δq, Q_q)

Return the precision-weighted squared source deviation.

# Arguments
- `Δq`, `Q_q`: source deviation and prior precision matrix

# Output
- `J_q`: scalar source-prior cost
"""
source_prior_cost(Δq, Q_q) = dot(Δq, Q_q, Δq)

"""
    mass_fraction_prior_cost(Δm, Q_m)

Return the precision-weighted squared mass-fraction deviation.

# Arguments
- `Δm`, `Q_m`: mass-fraction deviation and prior precision matrix

# Output
- `J_m`: scalar mass-fraction-prior cost
"""
mass_fraction_prior_cost(Δm, Q_m) = dot(Δm, Q_m, Δm)

"""
    costfunction(x, inversion)
    costfunction(x, inversion, workspace)
    costfunction(x, inversion, boundary_conditions, sources, mass_fractions)

Evaluate the complete inversion objective with allocated or cached controls.

# Arguments
- `x`, `inversion`: optimizer vector and inversion definition
- `workspace`, `boundary_conditions`, `sources`, `mass_fractions`: control values

# Output
- `J`: scalar sum of all objective components
"""
function costfunction(x::AbstractVector, inversion::Inversion)
    controls = unvec(inversion, x)
    return costfunction(x, inversion, controls.boundary_conditions,
        controls.sources, controls.mass_fractions)
end
function costfunction(x::AbstractVector, inversion::Inversion,
    workspace::NamedTuple)
    ranges = inversion.control_ranges
    workspace.boundary_values[inversion.boundary_control_indices] .=
        x[ranges.boundary]
    unvec!(workspace.boundary_conditions, workspace.boundary_values)
    unvec!(workspace.source_adjustments, x[ranges.source])
    unvec!(workspace.mass_fractions, x[ranges.mass_fraction],
        inversion.control_priors.mass_fraction)
    sources = adjustsource(inversion.control_priors.source,
        workspace.source_adjustments, inversion.control_templates.sources)
    return costfunction(x, inversion, workspace.boundary_conditions, sources,
        workspace.mass_fractions)
end
function costfunction(x::AbstractVector, inversion::Inversion,
    boundary_conditions::NamedTuple, sources::NamedTuple,
    mass_fractions::NamedTuple)
    observations = inversion.observations
    A = watermassmatrix(mass_fractions, observations.γ)
    ŷ = steadyinversion(lu(A), boundary_conditions, sources, observations;
        stoichiometry=inversion.stoichiometry)
    r = vec(ŷ) - vec(observations)
    Δx = x - inversion.x_prior
    ranges = inversion.control_ranges
    Δb = Δx[ranges.boundary]
    Δq = Δx[ranges.source]
    Δm = Δx[ranges.mass_fraction]
    J_y = data_cost(r, observations.Q)
    J_b = boundary_prior_cost(Δb, inversion.Q_b)
    J_s = boundary_smoothness_cost(Δb, inversion.Q_s)
    J_q = source_prior_cost(Δq, inversion.Q_q)
    J_m = mass_fraction_prior_cost(Δm, inversion.Q_m)
    return J_y + J_b + J_s + J_q + J_m
end

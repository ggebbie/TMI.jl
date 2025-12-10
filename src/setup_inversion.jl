include("./controls/control_parameters.jl")
include("./controls/boundary.jl")
include("./controls/massfrac.jl")
include("./controls/source.jl")
include("./controls/utils.jl")

"""
    setup_inversion(γ::Grid;
                    boundary::NamedTuple = NamedTuple(),
                    source::NamedTuple = NamedTuple(),
                    mass_fraction::NamedTuple = NamedTuple()
                   ) -> ControlParameters

A user-friendly function to set up the `ControlParameters` for an inversion problem.
This function abstracts away the complex internal structure of `ControlParameters`
and handles the creation of initial guesses, priors, and precision matrices.

# Arguments
- `γ::Grid`: The TMI Grid object.

# Keywords
Each keyword (`boundary`, `source`, `mass_fraction`) accepts a `NamedTuple` with the following optional fields:

- `prior`: (Required for the component to be active) A `NamedTuple` of `BoundaryCondition`, `Source`, or `Vector` objects representing the prior state.
- `initial_guess`: (Optional) The starting values for the control variables. Defaults to a `deepcopy` of the `prior`.
- `variance`: (Optional) A `NamedTuple` of scalar variances for each tracer/component. Used to construct diagonal precision matrices.
- `covariance`: (Optional) A `NamedTuple` of full covariance matrices. These will be inverted to form precision matrices.
- `lower_bound`, `upper_bound`: (Optional) `NamedTuple` of scalar bounds for each tracer/component.
- `couplings`: (Optional, for `source` only) A `NamedTuple` defining relationships between dependent and independent sources.

# Returns
- `controls::ControlParameters`: An initialized `ControlParameters` object ready for optimization.
```
"""
function setup_inversion(γ::Grid;
                         boundary::NamedTuple = NamedTuple(),
                         source::NamedTuple = NamedTuple(),
                         mass_fraction::NamedTuple = NamedTuple()
                        )

    # --- Boundary Conditions Setup ---
    u₀_b = get(boundary, :prior, NamedTuple())
    ub_ig = deepcopy(get(boundary, :initial_guess, u₀_b))
    ub_controls = !isempty(u₀_b) ? (isempty(ub_ig) ? deepcopy(u₀_b) : ub_ig) : NamedTuple()
    Qᵤ_final = !isempty(u₀_b) ? _build_tracer_precision_matrix(boundary, ub_controls) : NamedTuple()

    # --- Source Conditions Setup ---
    q₀_all = get(source, :prior, NamedTuple())
    couplings = get(source, :couplings, NamedTuple())
    
    uq_controls_nt = NamedTuple()
    if !isempty(q₀_all)
        dependent_sources = keys(couplings)
        
        uncertainty_config = get(source, :covariance, get(source, :variance, nothing))
        independent_controls_names = isnothing(uncertainty_config) ? () : keys(uncertainty_config)

        # Build the `uq` tuple, which has `nothing` for non-controlled sources
        potential_ind_sources = filter(k -> !in(k, dependent_sources), keys(q₀_all))
        uq_ig_full = deepcopy(get(source, :initial_guess, q₀_all))

        uq_controls_nt = NamedTuple{potential_ind_sources}(
             map(potential_ind_sources) do name
                if name in independent_controls_names
                    # Use initial guess if provided, otherwise default to the prior value.
                    get(uq_ig_full, name, q₀_all[name])
                else
                    # This source is not being controlled.
                    nothing
                end
            end
        )
    end
    
    uq_controls = isempty(uq_controls_nt) ? nothing : uq_controls_nt
    Qₛ_final = !isnothing(uq_controls) ? _build_tracer_precision_matrix(source, uq_controls) : NamedTuple()

    # --- Mass Fraction Setup ---
    m₀_mf = get(mass_fraction, :prior, NamedTuple())
    m_ig = deepcopy(get(mass_fraction, :initial_guess, m₀_mf))
    m_controls = !isempty(m₀_mf) ? (isempty(m_ig) ? deepcopy(m₀_mf) : m_ig) : NamedTuple()
    Qₘ_final = !isempty(m₀_mf) ? _build_massfrac_precision_matrix(mass_fraction, m_controls) : NamedTuple()

    # --- Bounds Setup ---
    raw_lb_b = get(boundary, :lower_bound, nothing)
    raw_ub_b = get(boundary, :upper_bound, nothing)
    lb_b_controls = _generate_control_bounds(ub_controls, raw_lb_b, -Inf)
    ub_b_controls = _generate_control_bounds(ub_controls, raw_ub_b, +Inf)

    raw_lb_q = get(source, :lower_bound, nothing)
    raw_ub_q = get(source, :upper_bound, nothing)
    lb_q_controls = !isnothing(uq_controls) ? _generate_control_bounds(uq_controls, raw_lb_q, -Inf) : NamedTuple()
    ub_q_controls = !isnothing(uq_controls) ? _generate_control_bounds(uq_controls, raw_ub_q, +Inf) : NamedTuple()

    raw_lb_m = get(mass_fraction, :lower_bound, nothing)
    raw_ub_m = get(mass_fraction, :upper_bound, nothing)
    if isnothing(raw_lb_m) && isnothing(raw_ub_m) && !isempty(m_controls)
        @warn "No explicit bounds provided for mass fractions. Defaulting to (-Inf, +Inf)."
    end
    lb_m_controls = _generate_control_bounds(m_controls, raw_lb_m, -Inf)
    ub_m_controls = _generate_control_bounds(m_controls, raw_ub_m, +Inf)
    
    lower_bound_vec = vectorize_controls(lb_b_controls, lb_q_controls, lb_m_controls)
    upper_bound_vec = vectorize_controls(ub_b_controls, ub_q_controls, ub_m_controls)

    return ControlParameters(; γ = γ,
                             ub = ub_controls, uq = uq_controls, m = m_controls,
                             u₀ = u₀_b, q₀ = q₀_all, m₀ = m₀_mf,
                             Qᵤ = Qᵤ_final, Qₛ = Qₛ_final, Qₘ = Qₘ_final,
                             lower_bound=lower_bound_vec, upper_bound=upper_bound_vec,
                             source_couplings = couplings)
end

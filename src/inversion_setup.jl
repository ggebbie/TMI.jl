
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

    # --- Set up the main control structures ---
    boundary_controls = BoundaryControls(boundary)
    source_controls = SourceControls(source)
    massfrac_controls = MassFracControls(mass_fraction, γ)

    # --- Bounds Setup ---
    raw_lb_b = get(boundary, :lower_bound, nothing)
    raw_ub_b = get(boundary, :upper_bound, nothing)
    lb_b_controls = _generate_control_bounds(boundary_controls.ub, raw_lb_b, -Inf)
    ub_b_controls = _generate_control_bounds(boundary_controls.ub, raw_ub_b, +Inf)

    raw_lb_q = get(source, :lower_bound, nothing)
    raw_ub_q = get(source, :upper_bound, nothing)
    lb_q_controls = _generate_control_bounds(source_controls.uq, raw_lb_q, -Inf)
    ub_q_controls = _generate_control_bounds(source_controls.uq, raw_ub_q, +Inf)

    raw_lb_m = get(mass_fraction, :lower_bound, nothing)
    raw_ub_m = get(mass_fraction, :upper_bound, nothing)
    if isnothing(raw_lb_m) && isnothing(raw_ub_m) && !isempty(massfrac_controls.m)
        @warn "No explicit bounds provided for mass fractions. Defaulting to (-Inf, +Inf)."
    end
    lb_m_controls = _generate_control_bounds(massfrac_controls.m, raw_lb_m, -Inf)
    ub_m_controls = _generate_control_bounds(massfrac_controls.m, raw_ub_m, +Inf)
    
    lower_bound_vec = vectorize_controls(lb_b_controls, lb_q_controls, lb_m_controls)
    upper_bound_vec = vectorize_controls(ub_b_controls, ub_q_controls, ub_m_controls)

    return ControlParameters(
        boundary_controls,
        source_controls,
        massfrac_controls,
        lower_bound_vec,
        upper_bound_vec,
        γ
    )
end

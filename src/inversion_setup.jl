
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
    boundary_controls = BoundaryControls(; boundary...)
    source_controls = SourceControls(; source...)
    massfrac_controls = MassFracControls(; mass_fraction..., γ=γ)

    # --- Vectorize Bounds ---
    lower_bound_vec = vectorize_controls(boundary_controls.lower_bound, source_controls.lower_bound, massfrac_controls.lower_bound)
    upper_bound_vec = vectorize_controls(boundary_controls.upper_bound, source_controls.upper_bound, massfrac_controls.upper_bound)

    return ControlParameters(
        boundary_controls,
        source_controls,
        massfrac_controls,
        lower_bound_vec,
        upper_bound_vec,
        γ
    )
end


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
Each keyword (`boundary`, `source`, `mass_fraction`) accepts a `NamedTuple`. If the tuple
is not empty, it MUST contain a `prior` field. Other optional fields include:

- `prior`: (Required) A `NamedTuple` of `BoundaryCondition`, `Source`, or `Vector` objects representing the prior state. This is mapped to `u₀`, `q₀`, or `m₀` internally as a positional argument.
- `initial_guess`: (Optional) The starting values for the control variables. This is mapped to `ub`, `uq`, or `m` internally. Defaults to a `deepcopy` of the `prior`.
- `variance`: (Optional) A `NamedTuple` of scalar variances for each tracer/component. Used to construct diagonal precision matrices.
- `covariance`: (Optional) A `NamedTuple` of full covariance matrices. These will be inverted to form precision matrices.
- `lower_bound`, `upper_bound`: (Optional) `NamedTuple` of scalar bounds for each tracer/component.
- `dependencies`: (Optional, for `source` only) A `NamedTuple` defining relationships between dependent and independent sources.

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
    if isempty(boundary)
        boundary_controls = BoundaryControls(nothing)
    else
        boundary_controls = BoundaryControls(
            boundary.prior,
            ub=get(boundary, :initial_guess, nothing),
            variance=get(boundary, :variance, nothing),
            covariance=get(boundary, :covariance, nothing),
            lower_bound=get(boundary, :lower_bound, nothing),
            upper_bound=get(boundary, :upper_bound, nothing)
        )
    end
    
    if isempty(source)
        source_controls = SourceControls(nothing)
    else
        source_controls = SourceControls(
            source.prior,
            uq=get(source, :initial_guess, nothing),
            variance=get(source, :variance, nothing),
            covariance=get(source, :covariance, nothing),
            dependencies=get(source, :dependencies, NamedTuple()),
            lower_bound=get(source, :lower_bound, nothing),
            upper_bound=get(source, :upper_bound, nothing)
        )
    end
    
    if isempty(mass_fraction)
        massfrac_controls = MassFracControls(nothing)
    else
        massfrac_controls = MassFracControls(
            mass_fraction.prior,
            m=get(mass_fraction, :initial_guess, nothing),
            variance=get(mass_fraction, :variance, nothing),
            covariance=get(mass_fraction, :covariance, nothing),
            lower_bound=get(mass_fraction, :lower_bound, nothing),
            upper_bound=get(mass_fraction, :upper_bound, nothing),
            γ=γ
        )
    end

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

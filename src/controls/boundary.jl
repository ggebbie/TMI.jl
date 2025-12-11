
"""
    struct BoundaryControls{U, U0, QU, D, G, B, LB, UB}

A container for control parameters related to model boundary conditions. This
includes the boundary conditions being optimized, their prior estimates, error
covariances, and cached values for gradients and perturbations.

# Fields
- `ub`: The boundary conditions to be optimized, typically a `NamedTuple` of
        `BoundaryCondition` objects.
- `u₀`: The prior (first-guess) estimate for the boundary conditions.
- `Qᵤ`: The inverse error covariance matrix for the boundary condition priors.
- `dub`: A cache for storing perturbations (`ub - u₀`) to the boundary
         conditions.
- `gdub`: A cache for storing the gradient of the cost function with respect to
          boundary condition perturbations.
- `b`: A cache for storing the full boundary condition field, initialized from
       the prior `u₀`.
- `lower_bound`: A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: A `NamedTuple` of upper bounds for each control variable.
"""
struct BoundaryControls{U, U0, QU, D, G, B, LB, UB}
    ub::U
    u₀::U0
    Qᵤ::QU
    dub::D
    gdub::G
    b::B
    lower_bound::LB
    upper_bound::UB
end

"""
    BoundaryControls(u₀::Union{NamedTuple, Nothing};
        ub=nothing,
        variance=nothing,
        covariance=nothing,
        lower_bound=nothing,
        upper_bound=nothing
    )

Constructs a `BoundaryControls` object. `u₀` is a required positional argument, but can be `nothing`.

If `u₀` is `nothing` or an empty `NamedTuple`, a null `BoundaryControls` object is returned where all fields are `nothing`.

# Arguments
- `u₀`: (Required) A `NamedTuple` of `BoundaryCondition` objects for the prior state, or `nothing`.
- `ub`: (Optional) The starting values for the control variables. Defaults to a `deepcopy` of `u₀`.
- `variance`: (Optional) A `NamedTuple` of scalar variances for each tracer.
- `covariance`: (Optional) A `NamedTuple` of full covariance matrices.
- `lower_bound`: (Optional) A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: (Optional) A `NamedTuple` of upper bounds for each control variable.
"""
function BoundaryControls(u₀::Union{NamedTuple, Nothing};
    ub=nothing,
    variance=nothing,
    covariance=nothing,
    lower_bound=nothing,
    upper_bound=nothing
)
    
    if isnothing(u₀) || isempty(u₀)
        @warn "No boundary controls (u₀) provided. Creating a null BoundaryControls object."
        return BoundaryControls(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end

    ub_controls = isnothing(ub) ? deepcopy(u₀) : deepcopy(ub)
    
    Qᵤ = _build_tracer_precision_matrix(ub_controls, variance, covariance)

    check_shared_references(ub_controls, "ub")
    check_shared_references(u₀, "u₀")
    check_shared_references(Qᵤ, "Qᵤ")

    # Bounds
    lower = _generate_control_bounds(ub_controls, lower_bound, -Inf)
    upper = _generate_control_bounds(ub_controls, upper_bound, +Inf)
    
    return BoundaryControls(
        ub_controls,
        u₀,
        Qᵤ,
        deepcopy(ub_controls), # dub
        deepcopy(ub_controls), # gdub
        deepcopy(u₀),  # b
        lower,
        upper
    )
end

"""
    update_b!(bc::BoundaryControls)

Update the main boundary condition buffer `bc.b` based on the current boundary
control variables.

This function applies the optimized boundary condition perturbations (`bc.ub`)
to the prior boundary conditions (`bc.u₀`) and stores the resulting boundary
conditions in `bc.b`. This `bc.b` then represents the full set of boundary
conditions used in the forward model run.
"""
function update_b!(bc::BoundaryControls)
    # Ensure bc.b is initialized with prior boundary conditions
    # Assuming bc.b has the same keys as bc.u₀, etc.

    du_names = keys(bc.ub) # Names of tracers for which boundary conditions are controlled
    @inbounds for key in keys(bc.b) # Iterate over all possible boundary condition keys
        setboundarycondition!(bc.b[key], bc.u₀[key]) # Start with the prior boundary condition
        if key ∈ du_names # If this tracer's boundary condition is being optimized
            # Calculate the perturbation: dub = ub - u₀
            setboundarycondition!(bc.dub[key], bc.ub[key])
            adjustboundarycondition!(bc.dub[key], bc.u₀[key]; r = -1.0)
            # Apply the perturbation to bc.b: bc.b = bc.b + dub
            adjustboundarycondition!(bc.b[key], bc.dub[key]; r = 1.0)
        end
    end
    return nothing
end

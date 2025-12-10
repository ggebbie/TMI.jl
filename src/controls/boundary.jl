
"""
    struct BoundaryControls{U, U0, QU, D, G, B}

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
"""
struct BoundaryControls{U, U0, QU, D, G, B}
    ub::U
    u₀::U0
    Qᵤ::QU
    dub::D
    gdub::G
    b::B
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

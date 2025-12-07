
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

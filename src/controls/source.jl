
"""
    struct SourceControls{Q, Q0, QS, D, G, S}

A container for control parameters related to interior sources and sinks. This
includes the source/sink terms being optimized, their prior estimates, error
covariances, and cached values for gradients and perturbations.

# Fields
- `uq`: The source/sink terms to be optimized, typically a `NamedTuple` of
        `Source` objects.
- `q₀`: The prior (first-guess) estimate for the source/sink terms.
- `Qₛ`: The inverse error covariance matrix for the source/sink priors.
- `duq`: A cache for storing perturbations (`uq - q₀`) to the source/sink terms.
- `gduq`: A cache for storing the gradient of the cost function with respect to
          source/sink perturbations.
- `q`: A cache for storing the full source/sink field, initialized from the
       prior `q₀`.
"""
struct SourceControls{Q, Q0, QS, D, G, S}
    uq::Q
    q₀::Q0
    Qₛ::QS
    duq::D
    gduq::G
    q::S
end

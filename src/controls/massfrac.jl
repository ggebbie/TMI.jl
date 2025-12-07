
"""
    struct MassFracControls{M, M0, QM, G, S, A}

A container for control parameters related to water-mass fractions, which
determine the transport matrix. This includes the mass fractions being optimized,
their prior estimates, error covariances, and cached values for gradients and
the transport matrix itself.

# Fields
- `m`: The mass fractions to be optimized, typically a `NamedTuple` of
       `MassFraction` objects.
- `m₀`: The prior (first-guess) estimate for the mass fractions.
- `Qₘ`: The inverse error covariance matrix for the mass fraction priors.
- `gm`: A cache for storing the gradient of the cost function with respect to
        the mass fractions.
- `steps`: A cache for precomputed step targets used in constructing the
           transport matrix.
- `A`: A cache for the assembled water-mass transport matrix.
"""
struct MassFracControls{M, M0, QM, G, S, A}
    m::M
    m₀::M0
    Qₘ::QM
    gm::G
    steps::S
    A::A
end

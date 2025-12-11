
"""
    struct MassFracControls{M, M0, QM, G, S, A, LB, UB}

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
- `lower_bound`: A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: A `NamedTuple` of upper bounds for each control variable.
"""
struct MassFracControls{M, M0, QM, G, S, A, LB, UB}
    m::M
    m₀::M0
    Qₘ::QM
    gm::G
    steps::S
    A::A
    lower_bound::LB
    upper_bound::UB
end

"""
    MassFracControls(;
        prior=NamedTuple(),
        initial_guess=nothing,
        variance=nothing,
        covariance=nothing,
        lower_bound=nothing,
        upper_bound=nothing,
        γ::Grid
    )

Constructs a `MassFracControls` object using keyword arguments.

# Arguments
- `prior`: (Required) A `NamedTuple` of `MassFraction` objects representing the prior state.
- `initial_guess`: (Optional) The starting values for the control variables. Defaults to a `deepcopy` of the `prior`.
- `variance`: (Optional) A `NamedTuple` of scalar variances for each component.
- `covariance`: (Optional) A `NamedTuple` of full covariance matrices.
- `lower_bound`: (Optional) A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: (Optional) A `NamedTuple` of upper bounds for each control variable.
- `γ`: The TMI grid object.
"""
function MassFracControls(;
    prior::NamedTuple=NamedTuple(),
    initial_guess=nothing,
    variance=nothing,
    covariance=nothing,
    lower_bound=nothing,
    upper_bound=nothing,
    γ::Grid
)
    if isempty(prior)
        return MassFracControls(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end

    m₀ = prior
    ig = isnothing(initial_guess) ? m₀ : initial_guess
    m = deepcopy(ig)
    
    Qₘ = _build_massfrac_precision_matrix(m, variance, covariance)

    check_shared_references(m, "m")
    check_shared_references(m₀, "m₀")

    # Precompute steps and cache the transport matrix
    m_steps = precompute_mass_fraction_steps(m, γ)
    A_cached = watermassmatrix(m, γ, m_steps)

    # Bounds
    lower = _generate_control_bounds(m, lower_bound, 0.0)
    upper = _generate_control_bounds(m, upper_bound, 1.0)
    
    return MassFracControls(
        m,
        m₀,
        Qₘ,
        deepcopy(m), # gm
        m_steps,
        A_cached,
        lower,
        upper
    )
end

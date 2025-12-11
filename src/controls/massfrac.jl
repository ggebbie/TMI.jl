
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
    MassFracControls(m₀::Union{NamedTuple, Nothing};
        m=nothing,
        variance=nothing,
        covariance=nothing,
        lower_bound=nothing,
        upper_bound=nothing,
        γ::Grid
    )

Constructs a `MassFracControls` object. `m₀` is a required positional argument, but can be `nothing`.

If `m₀` is `nothing` or an empty `NamedTuple`, a null `MassFracControls` object is returned where all fields are `nothing`.

# Arguments
- `m₀`: (Required) A `NamedTuple` of `MassFraction` objects representing the prior state, or `nothing`.
- `m`: (Optional) The starting values for the control variables. Defaults to a `deepcopy` of `m₀`.
- `variance`: (Optional) A scalar, `NamedTuple`, or vector of variances.
- `covariance`: (Optional) A full covariance matrix.
- `lower_bound`: (Optional) A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: (Optional) A `NamedTuple` of upper bounds for each control variable.
- `γ`: The TMI grid object.
"""
function MassFracControls(m₀::Union{NamedTuple, Nothing};
    m=nothing,
    variance=nothing,
    covariance=nothing,
    lower_bound=nothing,
    upper_bound=nothing,
    γ::Grid
)
    if isnothing(m₀) || isempty(m₀)
        @warn "No mass fraction controls (m₀) provided. Creating a null MassFracControls object."
        return MassFracControls(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end

    m_controls = isnothing(m) ? deepcopy(m₀) : deepcopy(m)
    
    Qₘ = _build_massfrac_precision_matrix(m_controls, variance, covariance)

    check_shared_references(m_controls, "m")
    check_shared_references(m₀, "m₀")

    # Precompute steps and cache the transport matrix
    m_steps = precompute_mass_fraction_steps(m_controls, γ)
    A_cached = watermassmatrix(m_controls, γ, m_steps)

    # Bounds
    lower = _generate_control_bounds(m_controls, lower_bound, 0.0)
    upper = _generate_control_bounds(m_controls, upper_bound, 1.0)
    
    return MassFracControls(
        m_controls,
        m₀,
        Qₘ,
        deepcopy(m_controls), # gm
        m_steps,
        A_cached,
        lower,
        upper
    )
end

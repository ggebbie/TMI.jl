
"""
    struct SourceControls{Q, Q0, QS, D, G, S, C, GC, LB, UB}

A container for control parameters related to interior sources and sinks.

This struct holds the source terms being optimized, their priors, error
covariances, and cached buffers. It also stores optional coupling information
to define dependencies between different tracer sources.

# Fields
- `uq`: `NamedTuple` of independent (optimizable) source terms.
- `q₀`: `NamedTuple` of the prior (first-guess) estimates for all sources.
- `Qₛ`: `NamedTuple` of inverse error covariance matrices for independent
        sources.
- `duq`: Cache for storing perturbations (`uq - q₀`) to independent sources.
- `gduq`: Cache for storing the final gradients with respect to the independent
          source controls.
- `q`: Cache for storing the full source field for all tracers, including
       the application of any couplings.
- `couplings`: `NamedTuple` defining dependencies between sources, e.g.,
               `(O2 = (:PO4, -170.0),)`.
- `gduq_cache`: A pre-allocated buffer to hold gradients for *all* sources
                before the chain rule is applied for couplings.
- `lower_bound`: A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: A `NamedTuple` of upper bounds for each control variable.
"""
struct SourceControls{Q, Q0, QS, D, G, S, C, GC, LB, UB}
    uq::Q
    q₀::Q0
    Qₛ::QS
    duq::D
    gduq::G
    q::S
    couplings::C
    gduq_cache::GC
    lower_bound::LB
    upper_bound::UB
end

"""
    SourceControls(;
        prior=NamedTuple(),
        initial_guess=nothing,
        variance=nothing,
        covariance=nothing,
        couplings=NamedTuple(),
        lower_bound=nothing,
        upper_bound=nothing
    )

Constructs a `SourceControls` object using keyword arguments.

# Arguments
- `prior`: (Required) A `NamedTuple` of `Source` objects representing the prior state for all sources.
- `initial_guess`: (Optional) The starting values for the control variables. Defaults to a `deepcopy` of the `prior`.
- `variance`: (Optional) A `NamedTuple` of scalar variances for each independent tracer.
- `covariance`: (Optional) A `NamedTuple` of full covariance matrices for independent tracers.
- `couplings`: (Optional) A `NamedTuple` defining relationships between dependent and independent sources.
- `lower_bound`: (Optional) A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: (Optional) A `NamedTuple` of upper bounds for each control variable.
"""
function SourceControls(;
    prior::NamedTuple=NamedTuple(),
    initial_guess=nothing,
    variance=nothing,
    covariance=nothing,
    couplings::NamedTuple=NamedTuple(),
    lower_bound=nothing,
    upper_bound=nothing
)
    if isempty(prior)
        return SourceControls(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end

    q₀_all = prior
    dependent_sources = keys(couplings)
    
    # Determine the names of sources that are explicitly controlled (have uncertainty)
    uncertainty_vars = if !isnothing(covariance)
        keys(covariance)
    elseif !isnothing(variance)
        keys(variance)
    else
        ()
    end

    # Independent sources are those that are not dependent on others
    potential_ind_sources = filter(k -> !in(k, dependent_sources), keys(q₀_all))
    
    # Use initial guess if provided, otherwise default to the prior value
    uq_ig_full = isnothing(initial_guess) ? q₀_all : initial_guess

    # Build the `uq` tuple, which has `nothing` for non-controlled sources
    uq_controls_nt = NamedTuple{potential_ind_sources}(
        map(potential_ind_sources) do name
            if name in uncertainty_vars
                # This is an independent control variable
                deepcopy(get(uq_ig_full, name, q₀_all[name]))
            else
                # This source is not being controlled
                nothing
            end
        end
    )
    
    uq_controls = isempty(uq_controls_nt) ? nothing : uq_controls_nt
    Qₛ_final = !isnothing(uq_controls) ? _build_tracer_precision_matrix(uq_controls, variance, covariance) : NamedTuple()

    independent_sources = isnothing(uq_controls) ? () : filter(k -> !isnothing(uq_controls[k]), keys(uq_controls))
    uq_indep = isnothing(uq_controls) ? nothing : NamedTuple{independent_sources}(getproperty(uq_controls, k) for k in independent_sources)
    Qₛ_indep = isnothing(Qₛ_final) ? nothing : NamedTuple{independent_sources}(getproperty(Qₛ_final, k) for k in independent_sources)

    check_shared_references(uq_indep, "uq")
    check_shared_references(q₀_all, "q₀")
    check_shared_references(Qₛ_indep, "Qₛ")

    if !isnothing(uq_indep)
        for dep_name in keys(couplings)
            if haskey(uq_indep, dep_name) && !isnothing(uq_indep[dep_name])
                error("Tracer :$(dep_name) is dependent in `source_couplings`\n" *
                      "and independent in `uq`. Source cannot be both.")
            end
        end
    end

    # Bounds
    lower = _generate_control_bounds(uq_indep, lower_bound, -Inf)
    upper = _generate_control_bounds(uq_indep, upper_bound, +Inf)

    return SourceControls(
        uq_indep,
        q₀_all,
        Qₛ_indep,
        isnothing(uq_indep) ? nothing : deepcopy(uq_indep), # duq
        isnothing(uq_indep) ? nothing : deepcopy(uq_indep), # gduq
        deepcopy(q₀_all), # q
        couplings,
        zero(q₀_all), # gduq_cache
        lower,
        upper
    )
end


"""
    update_q!(sc::SourceControls)

Update the main source buffer `sc.q` and perturbation buffer `sc.duq`.

This function orchestrates the application of source controls by:
1. Initializing the buffer `sc.q` with the prior sources `sc.q₀`.
2. Applying the independent (optimized) source perturbations from `sc.uq`.
3. Calculating and applying any coupled sources defined in `sc.couplings`.
"""
function update_q!(sc::SourceControls)
    # Return early if there are no sources to update
    isnothing(sc.q) && return nothing

    # 1. Initialize q from the prior q₀
    replacesource!(sc.q, sc.q₀)

    # 2. Apply independent source controls from uq
    if !isnothing(sc.uq)
        for key in keys(sc.uq)
            # Update perturbation duq = uq - q₀
            replacesource!(sc.duq[key], sc.uq[key])
            adjustsource!(sc.duq[key], sc.q₀[key]; r=-1.0)

            # Update main buffer q = q₀ + duq
            adjustsource!(sc.q[key], sc.duq[key])
        end
    end

    # 3. Apply coupled source logic
    if !isnothing(sc.couplings)
        for dep_name in keys(sc.couplings)
            indep_name, alpha = sc.couplings[dep_name]
            
            # Get the final source of the independent tracer
            independent_source_field = sc.q[indep_name]

            # Calculate and overwrite the dependent source
            coupled_source_field = alpha * independent_source_field
            replacesource!(sc.q[dep_name], coupled_source_field)
        end
    end

    return nothing
end

"""
    update_gduq!(sc::SourceControls)

Adjust the source gradient buffer `sc.gduq` to account for all gradient
contributions, including priors, direct misfit, and coupled sources.

This function assumes `sc.gduq` already contains the prior gradient. It adds
the direct model-data misfit gradient from `sc.gduq_cache`. Then, if
couplings are defined, it applies the chain rule to accumulate the scaled
gradients from any dependent tracers.
"""
function update_gduq!(sc::SourceControls)
    # sc.gduq already contains the prior gradients.
    # sc.gduq_cache contains the misfit gradients for ALL sources.

    # 1. Add the direct misfit gradients to the independent controls' gradients.
    #    This performs: gduq = gduq (from prior) + gduq_cache (from misfit)
    for key in keys(sc.gduq)
        adjustsource!(sc.gduq[key], sc.gduq_cache[key])
    end

    # 2. If couplings exist, add the coupled misfit gradients via the chain rule.
    if !isnothing(sc.couplings)
        for dep_name in keys(sc.couplings)
            indep_name, alpha = sc.couplings[dep_name]
            
            if haskey(sc.gduq, indep_name)
                 # Accumulate: gduq[i] += alpha * gduq_cache[j]
                 adjustsource!(sc.gduq[indep_name], sc.gduq_cache[dep_name]; r=alpha)
            end
        end
    end

    return nothing
end

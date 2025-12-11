
"""
    struct SourceControls{Q, Q0, QS, D, G, S, C, GC, LB, UB}

A container for control parameters related to interior sources and sinks.

This struct holds the source terms being optimized, their priors, error
covariances, and cached buffers. It also stores optional dependency information
to define relationships between different tracer sources.

# Fields
- `uq`: `NamedTuple` of independent (optimizable) source terms.
- `q₀`: `NamedTuple` of the prior (first-guess) estimates for all sources.
- `Qₛ`: `NamedTuple` of inverse error covariance matrices for independent
        sources.
- `duq`: Cache for storing perturbations (`uq - q₀`) to independent sources.
- `gduq`: Cache for storing the final gradients with respect to the independent
          source controls.
- `q`: Cache for storing the full source field for all tracers, including
       the application of any dependencies.
- `dependencies`: `NamedTuple` defining dependencies between sources, e.g.,
                  `(O2 = (:PO4, -170.0),)`.
- `gduq_cache`: A pre-allocated buffer to hold gradients for *all* sources
                before the chain rule is applied for dependencies.
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
    dependencies::C
    gduq_cache::GC
    lower_bound::LB
    upper_bound::UB
end

"""
    SourceControls(q₀::Union{NamedTuple, Nothing};
        uq=nothing,
        variance=nothing,
        covariance=nothing,
        dependencies=NamedTuple(),
        lower_bound=nothing,
        upper_bound=nothing
    )

Constructs a `SourceControls` object. `q₀` is a required positional argument, but can be `nothing`.

If `q₀` is `nothing` or an empty `NamedTuple`, a null `SourceControls` object is returned where all fields are `nothing`.

# Arguments
- `q₀`: (Required) A `NamedTuple` of `Source` objects representing the prior state for all sources, or `nothing`.
- `uq`: (Optional) The starting values for the control variables (independent sources). Defaults to a `deepcopy` of `q₀`.
- `variance`: (Optional) A `NamedTuple` of scalar variances for each independent tracer.
- `covariance`: (Optional) A `NamedTuple` of full covariance matrices for independent tracers.
- `dependencies`: (Optional) A `NamedTuple` defining relationships between dependent and independent sources, of the form `(dependent = (:independent, scaling_factor),)`.
- `lower_bound`: (Optional) A `NamedTuple` of lower bounds for each control variable.
- `upper_bound`: (Optional) A `NamedTuple` of upper bounds for each control variable.
"""
function SourceControls(q₀::Union{NamedTuple, Nothing};
    uq=nothing,
    variance=nothing,
    covariance=nothing,
    dependencies::NamedTuple=NamedTuple(),
    lower_bound=nothing,
    upper_bound=nothing
)
    if isnothing(q₀) || isempty(q₀)
        @warn "No source controls (q₀) provided. Creating a null SourceControls object."
        return SourceControls(nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing, nothing)
    end

    dependent_sources = keys(dependencies)
    
    uncertainty_vars = if !isnothing(covariance)
        keys(covariance)
    elseif !isnothing(variance)
        keys(variance)
    else
        ()
    end

    potential_ind_sources = filter(k -> !in(k, dependent_sources), keys(q₀))
    
    uq_ig_full = isnothing(uq) ? q₀ : uq

    uq_controls_nt = NamedTuple{potential_ind_sources}(
        map(potential_ind_sources) do name
            if name in uncertainty_vars
                deepcopy(get(uq_ig_full, name, q₀[name]))
            else
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
    check_shared_references(q₀, "q₀")
    check_shared_references(Qₛ_indep, "Qₛ")

    if !isnothing(uq_indep)
        for dep_name in keys(dependencies)
            if haskey(uq_indep, dep_name) && !isnothing(uq_indep[dep_name])
                error("Tracer :$(dep_name) is listed in `dependencies`\n" *
                      "but is also an independent control variable in `uq`. A source cannot be both.")
            end
        end
    end

    lower = _generate_control_bounds(uq_indep, lower_bound, -Inf)
    upper = _generate_control_bounds(uq_indep, upper_bound, +Inf)

    return SourceControls(
        uq_indep,
        q₀,
        Qₛ_indep,
        isnothing(uq_indep) ? nothing : deepcopy(uq_indep),
        isnothing(uq_indep) ? nothing : deepcopy(uq_indep),
        deepcopy(q₀),
        dependencies,
        zero(q₀),
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
3. Calculating and applying any dependent sources defined in `sc.dependencies`.
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

    # 3. Apply dependent source logic
    if !isnothing(sc.dependencies)
        for dep_name in keys(sc.dependencies)
            indep_name, alpha = sc.dependencies[dep_name]
            
            # Get the final source of the independent tracer
            independent_source_field = sc.q[indep_name]

            # Calculate and overwrite the dependent source
            dependent_source_field = alpha * independent_source_field
            replacesource!(sc.q[dep_name], dependent_source_field)
        end
    end

    return nothing
end

"""
    update_gduq!(sc::SourceControls)

Adjust the source gradient buffer `sc.gduq` to account for all gradient
contributions, including priors, direct misfit, and dependent sources.

This function assumes `sc.gduq` already contains the prior gradient. It adds
the direct model-data misfit gradient from `sc.gduq_cache`. Then, if
dependencies are defined, it applies the chain rule to accumulate the scaled
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

    # 2. If dependencies exist, add the dependent misfit gradients via the chain rule.
    if !isnothing(sc.dependencies)
        for dep_name in keys(sc.dependencies)
            indep_name, alpha = sc.dependencies[dep_name]
            
            if haskey(sc.gduq, indep_name)
                 # Accumulate: gduq[i] += alpha * gduq_cache[j]
                 adjustsource!(sc.gduq[indep_name], sc.gduq_cache[dep_name]; r=alpha)
            end
        end
    end

    return nothing
end

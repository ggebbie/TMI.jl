
"""
    struct SourceControls{Q, Q0, QS, D, G, S, C, GC}

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
"""
struct SourceControls{Q, Q0, QS, D, G, S, C, GC}
    uq::Q
    q₀::Q0
    Qₛ::QS
    duq::D
    gduq::G
    q::S
    couplings::C
    gduq_cache::GC
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

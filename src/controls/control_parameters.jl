

"""
    struct ControlParameters{BC, SC, MC, V, LB, UB, UBL, QL, ML, G <: Grid}

A top-level container for all control variables in an optimization problem. It
aggregates boundary, source, and mass fraction controls, along with the grid,
cached data, and optimization bounds for efficient computation.
"""
struct ControlParameters{BC, SC, MC, V, LB, UB, UBL, QL, ML, G <: Grid}
    boundary::BC
    source::SC
    massfrac::MC
    vector::V
    lower_bound::LB
    upper_bound::UB
    ub_length::UBL
    uq_length::QL
    m_length::ML
    γ::G
end


"""
    ControlParameters(boundary_controls, source_controls, massfrac_controls, lower_bound, upper_bound, γ; cache_vector=false, cache_lengths=true)

Construct and initialize the `ControlParameters` object for an optimization.
This function takes all individual control components (boundary conditions,
sources, etc.), organizes them into the appropriate sub-structs, and performs
necessary pre-caching.
"""
function ControlParameters(boundary_controls::BoundaryControls,
                           source_controls::SourceControls,
                           massfrac_controls::MassFracControls,
                           lower_bound::AbstractVector,
                           upper_bound::AbstractVector,
                           γ::Grid;
                           cache_vector=false, cache_lengths=true)

    # --- Final Assembly ---
    ub_length = cache_lengths && !isnothing(boundary_controls.ub) ? length(vec(boundary_controls.ub)) : nothing
    uq_length = cache_lengths && !isnothing(source_controls.uq) ? length(vec(source_controls.uq)) : nothing
    m_length = cache_lengths && !isnothing(massfrac_controls.m) ? length(vec(massfrac_controls.m)) : nothing

    vector = cache_vector ? vectorize_controls(boundary_controls.ub,
                                               source_controls.uq,
                                               massfrac_controls.m) : nothing
    
    return ControlParameters(boundary_controls, source_controls, massfrac_controls,
                             vector, lower_bound, upper_bound,
                             ub_length, uq_length, m_length, γ)
end


"""
    vectorize_controls(ub, uq, m) -> Vector

Flatten the core control variables (`ub`, `uq`, `m`) into a single state
vector. This is a necessary step for interfacing with optimization algorithms
that expect a single vector of parameters.
"""
function vectorize_controls(ub, uq, m)
    v = []
    !isnothing(ub) && append!(v, vec(ub))
    !isnothing(uq) && append!(v, vec(uq))
    !isnothing(m) && append!(v, vec(m))
    return v
end

"""
    vec(cp::ControlParameters) -> Vector

Extract or compute the flattened vector representation of the control
parameters. If the vector has been cached in `cp.vector`, it is returned
directly; otherwise, it is computed by calling `vectorize_controls`.
"""
function vec(cp::ControlParameters)
    if !isnothing(cp.vector)
        return cp.vector
    else
        return vectorize_controls(cp.boundary.ub, cp.source.uq, cp.massfrac.m)
    end
end

"""
    unvec(cp::ControlParameters, vector::AbstractVector) -> (ub, uq, m)

Reconstruct the structured control variables (`ub`, `uq`, `m`) from a flat
state vector. This function uses the cached lengths in `cp` to correctly
partition the vector and reshape the components back into their original
`NamedTuple` format.
"""
function unvec(cp::ControlParameters, vector::AbstractVector)
    # Calculate boundaries
    u_len = cp.ub_length
    q_len = cp.uq_length
    m_len = cp.m_length
    
    idx = 1
    ub = nothing; uq = nothing; m = nothing
    
    if !isnothing(cp.boundary.ub)
        ub = unvec(cp.boundary.ub, vector[idx:idx+u_len-1])
        idx += u_len
    end
    
    if !isnothing(cp.source.uq)
        uq = unvec(cp.source.uq, vector[idx:idx+q_len-1])
        idx += q_len
    end
    
    if !isnothing(cp.massfrac.m)
        m = unvec(cp.massfrac.m, vector[idx:idx+m_len-1])
    end

    return ub, uq, m
end

"""
    unvec!(cp::ControlParameters, vector::AbstractVector)

Perform an in-place update of the control variables stored within `cp` from a
flat state vector. This is a performance-optimized version of `unvec` that
modifies the `ub`, `uq`, and `m` fields directly, avoiding extra memory
allocations.
"""
function unvec!(cp::ControlParameters, vector::AbstractVector)
    idx = 1

    if !isnothing(cp.boundary.ub)
        idx = unvec!(cp.boundary.ub, vector; idx=idx, return_idx=true)
    end

    if !isnothing(cp.source.uq)
        idx = unvec!(cp.source.uq, vector; idx=idx, return_idx=true)
    end

    if !isnothing(cp.massfrac.m)
        unvec!(cp.massfrac.m, vector; idx=idx)
    end

    return cp.boundary.ub, cp.source.uq, cp.massfrac.m
end




"""
    struct Controls{BC, SC, MC, V, LB, UB, UBL, QL, ML, G <: Grid}

A top-level container for all control variables in an optimization problem. It
aggregates boundary, source, and mass fraction controls, along with the grid,
cached data, and optimization bounds for efficient computation.
"""
struct Controls{BC, SC, MC, V, LB, UB, UBL, QL, ML, G <: Grid}
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
    Controls(γ::Grid;
        boundary::BoundaryControls = BoundaryControls(nothing),
        source::SourceControls = SourceControls(nothing),
        massfrac::MassFracControls = MassFracControls(nothing, γ=γ),
        cache_vector=false,
        cache_lengths=true
    ) -> Controls

A user-friendly, outer constructor that assembles the main `Controls` object from its constituent control-variable parts.

This constructor takes pre-built `BoundaryControls`, `SourceControls`, and `MassFracControls` objects, automatically extracts and vectorizes their bounds, and then calls the inner constructor to build the final `Controls` object.

# Arguments
- `γ::Grid`: The TMI Grid object.
- `boundary`: A `BoundaryControls` object.
- `source`: A `SourceControls` object.
- `massfrac`: A `MassFracControls` object.
- `cache_vector`, `cache_lengths`: (Optional) Flags for caching internal vector representations.

# Returns
- `controls::Controls`: An initialized `Controls` object ready for optimization.
"""
function Controls(γ::Grid;
    boundary::BoundaryControls = BoundaryControls(nothing),
    source::SourceControls = SourceControls(nothing),
    massfrac::MassFracControls = MassFracControls(nothing, γ=γ),
    cache_vector=false,
    cache_lengths=true
)
    # --- Vectorize Bounds ---
    lower_bound_vec = vectorize_controls(boundary.lower_bound, source.lower_bound, massfrac.lower_bound)
    upper_bound_vec = vectorize_controls(boundary.upper_bound, source.upper_bound, massfrac.upper_bound)

    return Controls(
        boundary,
        source,
        massfrac,
        lower_bound_vec,
        upper_bound_vec,
        γ,
        cache_vector=cache_vector,
        cache_lengths=cache_lengths
    )
end

"""
    Controls(boundary_controls, source_controls, massfrac_controls, lower_bound, upper_bound, γ; cache_vector=false, cache_lengths=true)

Construct and initialize the `Controls` object for an optimization.
This function takes all individual control components (boundary conditions,
sources, etc.), organizes them into the appropriate sub-structs, and performs
necessary pre-caching.
"""
function Controls(boundary_controls::BoundaryControls,
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
    
    return Controls(boundary_controls, source_controls, massfrac_controls,
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
    vec(cp::Controls) -> Vector

Extract or compute the flattened vector representation of the control
parameters. If the vector has been cached in `cp.vector`, it is returned
directly; otherwise, it is computed by calling `vectorize_controls`.
"""
function vec(cp::Controls)
    if !isnothing(cp.vector)
        return cp.vector
    else
        return vectorize_controls(cp.boundary.ub, cp.source.uq, cp.massfrac.m)
    end
end

"""
    unvec(cp::Controls, vector::AbstractVector) -> (ub, uq, m)

Reconstruct the structured control variables (`ub`, `uq`, `m`) from a flat
state vector. This function uses the cached lengths in `cp` to correctly
partition the vector and reshape the components back into their original
`NamedTuple` format.
"""
function unvec(cp::Controls, vector::AbstractVector)
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
    unvec!(cp::Controls, vector::AbstractVector)

Perform an in-place update of the control variables stored within `cp` from a
flat state vector. This is a performance-optimized version of `unvec` that
modifies the `ub`, `uq`, and `m` fields directly, avoiding extra memory
allocations.
"""
function unvec!(cp::Controls, vector::AbstractVector)
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

"""
    write_gradient!(G, gu, gq, gm, controls)
    Flatten structured gradients into a single vector.
# Arguments
- `G::Vector`: flat gradient vector to be populated in-place.
- `gu`: gradient of the boundary controls.
- `gq`: gradient of the source controls.
- `gm`: gradient of the mass-fraction controls.
- `controls::Controls`: container describing which control blocks are active and their shapes.
# Output
- `nothing`: operates in-place, mutating `G`.
"""
@inline function write_gradient!(G, gu, gq, gm, controls::Controls)
    idx = 1
    if !isnothing(controls.boundary.ub)
        guv = vec(gu)
        len = length(guv)
        copyto!(G, idx, guv, 1, len)
        idx += len
    end
    if !isnothing(controls.source.uq)
        gqv = vec(gq)
        len = length(gqv)
        copyto!(G, idx, gqv, 1, len)
        idx += len
    end
    if !isnothing(controls.massfrac.m)
        gmv = vec(gm)
        len = length(gmv)
        copyto!(G, idx, gmv, 1, len)
    end
    return nothing
end

"""
    descale_parameter(x::NamedTuple, scale::NamedTuple)

Return a deepcopy of `x` with each entry divided by its corresponding scale.
Supports numeric entries and `Field`/`BoundaryCondition`/`Source` via their
`tracer` arrays; `nothing` entries are left untouched.
"""
function descale_parameter(x::NamedTuple, scale::NamedTuple)
    x_copy = deepcopy(x)
    for k in keys(x)
        v = x[k]
        if !isnothing(v)
            if v isa Number
                x_copy = merge(x_copy, NamedTuple{(k,)}((v / scale[k],)))
            elseif v isa Union{Field, BoundaryCondition, Source} 
                x_copy[k].tracer ./= scale[k]
            else
                @error "X must only contain numbers" value=v key=k
            end
        end
    end
    return x_copy
end

"""
    descale_parameter!(x::NamedTuple, scale::NamedTuple)

In-place version of `descale_parameter`, dividing supported entries of `x`
by the provided per-key scale factors.
"""
function descale_parameter!(x::NamedTuple, scale::NamedTuple)
    for name in eachindex(x)
        if !isnothing(x[name])
            if x[name] isa Union{Field, BoundaryCondition, Source}
                x[name].tracer ./= scale[name]
            end
        end
    end
end

"""
    rescale_parameter!(x::NamedTuple, scale::NamedTuple)

In-place inverse of `descale_parameter!`: multiply supported entries of `x`
by the provided per-key scale factors.
"""
function rescale_parameter!(x::NamedTuple, scale::NamedTuple)
    for name in eachindex(x)
        if !isnothing(x[name])
            if x[name] isa Union{Field, BoundaryCondition} 
                x[name].tracer .*= scale[name]
            end
        end
    end
end

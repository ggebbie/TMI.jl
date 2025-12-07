
include("controls/boundary.jl")
include("controls/source.jl")
include("controls/massfrac.jl")

"""
    check_shared_references(nt::NamedTuple, name::String)

Check if any fields in a `NamedTuple` point to the same object in memory. This
is a utility to prevent bugs where in-place updates to one control parameter
accidentally modify another. Throws an error if shared references are found.
"""
function check_shared_references(nt::NamedTuple, name::String)
    vals = filter(!isnothing, collect(values(nt)))
    n = length(vals)

    for i in 1:n
        for j in (i+1):n
            if vals[i] === vals[j]
                key_names = collect(keys(nt))
                # Find which keys share the reference
                shared_keys = String[]
                for k in 1:length(nt)
                    if !isnothing(nt[k]) && (nt[k] === vals[i])
                        push!(shared_keys, string(key_names[k]))
                    end
                end
                error("Shared reference detected in $name: fields " *
                      "$(join(shared_keys, ", ")) point to the same object in " *
                      "memory. Use deepcopy to create independent copies.")
            end
        end
    end
end
check_shared_references(x, name::String) = nothing  # Non-NamedTuple, no check needed

"""
    struct ControlParameters{BC, SC, MC, V, UBL, QL, ML, G <: Grid}

A top-level container for all control variables in an optimization problem. It
aggregates boundary, source, and mass fraction controls, along with the grid
and cached data for efficient computation. This struct is the primary object
passed to cost functions and optimization routines.
"""
struct ControlParameters{BC, SC, MC, V, UBL, QL, ML, G <: Grid}
    boundary::BC
    source::SC
    massfrac::MC
    vector::V
    ub_length::UBL
    uq_length::QL
    m_length::ML
    γ::G
end

"""
    ControlParameters(; γ::Grid, ub=nothing, uq=nothing, m=nothing, ...)

Construct and initialize the `ControlParameters` object for an optimization.
This function takes all individual control components (boundary conditions,
sources, etc.), organizes them into the appropriate sub-structs
(`BoundaryControls`, `SourceControls`, `MassFracControls`), and performs
necessary pre-caching. Caching can include the transport matrix `A` and a
flattened vector of all control variables for use with optimizers.
"""
function ControlParameters(; γ::Grid,
                            ub=nothing, uq=nothing, m=nothing,
                            u₀=nothing, q₀=nothing, m₀=nothing,
                            Qᵤ=nothing, Qₛ=nothing, Qₘ=nothing,
                            cache_vector=false, cache_lengths=true,
                            cache_mass_fraction_steps=true,
                            cache_A=true)

    for (val, name) in ((ub, "ub"), (uq, "uq"), (m, "m"),
                        (u₀, "u₀"), (q₀, "q₀"), (m₀, "m₀"),
                        (Qᵤ, "Qᵤ"), (Qₛ, "Qₛ"))
        check_shared_references(val, name)
    end

    boundary_controls = BoundaryControls(
        ub, u₀, Qᵤ,
        isnothing(ub) ? nothing : deepcopy(ub),
        isnothing(ub) ? nothing : deepcopy(ub),
        isnothing(u₀) ? nothing : deepcopy(u₀)
    )

    source_controls = SourceControls(
        uq, q₀, Qₛ,
        isnothing(uq) ? nothing : deepcopy(uq),
        isnothing(uq) ? nothing : deepcopy(uq),
        isnothing(q₀) ? nothing : deepcopy(q₀)
    )
    
    m_steps = nothing
    if cache_mass_fraction_steps
        isnothing(m) && throw(ArgumentError("cache_mass_fraction_steps=true requires m"))
        m_steps = precompute_mass_fraction_steps(m, γ)
    end

    A_cached = nothing
    if cache_A
        isnothing(m) && throw(ArgumentError("cache_A=true requires m"))
        A_cached = watermassmatrix(m, γ, cache_mass_fraction_steps ? m_steps : nothing)
    end

    massfrac_controls = MassFracControls(
        m, m₀, Qₘ,
        isnothing(m) ? nothing : deepcopy(m),
        m_steps,
        A_cached
    )

    ub_length = cache_lengths && !isnothing(ub) ? length(vec(ub)) : nothing
    uq_length = cache_lengths && !isnothing(uq) ? length(vec(uq)) : nothing
    m_length = cache_lengths && !isnothing(m) ? length(vec(m)) : nothing

    vector = cache_vector ? vectorize_controls(boundary_controls.ub,
                                               source_controls.uq,
                                               massfrac_controls.m) : nothing
    
    return ControlParameters(boundary_controls, source_controls, massfrac_controls,
                             vector, ub_length, uq_length, m_length, γ)
end


"""
    vectorize_controls(ub, uq, m) -> Vector

Flatten the core control variables (`ub`, `uq`, `m`) into a single state
vector. This is a necessary step for interfacing with optimization algorithms
that expect a single vector of parameters.
"""
function vectorize_controls(ub::NamedTuple, uq::NamedTuple, m::NamedTuple)
    return vcat(vec.([ub, uq, m])...)
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

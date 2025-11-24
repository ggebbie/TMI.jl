"""
    ControlParameters{u_names, q_names, m_names, u₀_names, q₀_names, m₀_names}

Container for optimization control variables (boundary conditions, sources, mass fractions) with vectorization support.
Includes optional first guess values and weighting matrices for inverse problems.

# Fields
- `du`: NamedTuple of boundary condition perturbations (BoundaryCondition objects)
- `dq`: NamedTuple of source perturbations (Source objects)
- `m`: NamedTuple of mass fractions (MassFraction objects)
- `u₀`: NamedTuple of first guess boundary conditions
- `q₀`: NamedTuple of first guess sources
- `m₀`: NamedTuple of first guess mass fractions
- `Qⁱᵤ`: NamedTuple of inverse covariance matrices (Symmetric) for boundary condition errors
- `Qⁱₛ`: NamedTuple of inverse covariance matrices (Symmetric) for source errors
- `Qⁱ_m`: Single inverse covariance matrix (Symmetric) for mass fraction errors
- `vector`: Flattened vector representation of all control parameters
- `du_lengths`: Dictionary mapping field names to vector lengths for `du` (if store_lengths=true)
- `dq_lengths`: Dictionary mapping field names to vector lengths for `dq` (if store_lengths=true)
- `m_lengths`: Dictionary mapping field names to vector lengths for `m` (if store_lengths=true)

# Constructor
```julia
ControlParameters(; du = nothing, dq = nothing, m = nothing,
                   u₀ = nothing, q₀ = nothing, m₀ = nothing,
                   Qⁱᵤ = nothing, Qⁱₛ = nothing, Qⁱ_m = nothing,
                   store_lengths = false)
```

All parameters are optional and default to `nothing`. Set `store_lengths=true` to enable
unraveling vectors back to structured NamedTuples using the `unravel` function.

# Examples
```julia
# Minimal usage with only boundary conditions
cp = ControlParameters(du = (θ = bc1, S = bc2))

# With first guesses and weighting matrices
cp = ControlParameters(
    du = (θ = bc_pert,),
    u₀ = (θ = bc_guess,),
    Qⁱᵤ = (θ = inv_cov_matrix,),
    store_lengths = true
)

# Unravel control vector back to structured form
du_new, dq_new, m_new = unravel(cp, optimized_vector)
```

See also: [`unravel`](@ref), [`unravel_field`](@ref)
"""
struct ControlParameters{u_names, q_names, m_names, u₀_names, q₀_names, m₀_names}
    du::Union{NamedTuple{u_names, <:Tuple{Vararg{BoundaryCondition}}}, Nothing}
    dq::Union{NamedTuple{q_names, <:Tuple{Vararg{Source}}}, Nothing}
    m::Union{NamedTuple{m_names, <:Tuple{Vararg{MassFraction}}}, Nothing}
    u₀::Union{NamedTuple{u₀_names, <:Tuple{Vararg{BoundaryCondition}}}, Nothing}
    q₀::Union{NamedTuple{q₀_names, <:Tuple{Vararg{Source}}}, Nothing}
    m₀::Union{NamedTuple{m₀_names, <:Tuple{Vararg{MassFraction}}}, Nothing}
    Qⁱᵤ::Union{NamedTuple{u₀_names, <:Tuple{Vararg{Symmetric}}}, Nothing}
    Qⁱₛ::Union{NamedTuple{q₀_names, <:Tuple{Vararg{Symmetric}}}, Nothing}
    Qⁱ_m::Union{Symmetric, Nothing}
    vector::Vector
    du_lengths::Union{Dict{Symbol, Int}, Nothing}
    dq_lengths::Union{Dict{Symbol, Int}, Nothing}
    m_lengths::Union{Dict{Symbol, Int}, Nothing}

    function ControlParameters(; du = nothing, dq = nothing, m = nothing,
                                u₀ = nothing, q₀ = nothing, m₀ = nothing,
                                Qⁱᵤ = nothing, Qⁱₛ = nothing, Qⁱ_m = nothing,
                                store_lengths = false)
        vector = vcat(vec.(filter(!isnothing, [du, dq, m]))...)

        u_names = isnothing(du) ? () : typeof(du).parameters[1]
        q_names = isnothing(dq) ? () : typeof(dq).parameters[1]
        m_names = isnothing(m) ? () : typeof(m).parameters[1]
        u₀_names = isnothing(u₀) ? () : typeof(u₀).parameters[1]
        q₀_names = isnothing(q₀) ? () : typeof(q₀).parameters[1]
        m₀_names = isnothing(m₀) ? () : typeof(m₀).parameters[1]

        if store_lengths
            du_lengths = isnothing(du) ? 0 : Dict(k => length(vec(v)) for (k,v) in pairs(du))
            dq_lengths = isnothing(dq) ? 0 : Dict(k => length(vec(v)) for (k,v) in pairs(dq))
            m_lengths = isnothing(m) ? 0 : Dict(k => length(vec(v)) for (k,v) in pairs(m))
        else
            du_lengths = nothing
            dq_lengths = nothing
            m_lengths = nothing
        end

        new{u_names, q_names, m_names, u₀_names, q₀_names, m₀_names}(du, dq, m, u₀, q₀, m₀, Qⁱᵤ, Qⁱₛ, Qⁱ_m, vector, du_lengths, dq_lengths, m_lengths)
    end
end

"""
    unravel_field(field, lengths, vector) -> NamedTuple

Reconstruct a NamedTuple of typed fields from a flat vector using stored lengths.
"""
function unravel_field(
    field::NamedTuple{names, <:Tuple{Vararg{T}}},
    lengths::Union{Dict{Symbol,Int}, Nothing},
    vector::Vector) where {names, T}
    
    isnothing(lengths) && error("Lengths not stored. Create ControlParameters with store_lengths=true")
    
    n_fields = length(names)
    vals = Vector{T}(undef, n_fields)
    
    idx = 1
    for (i, k) in enumerate(names)
        n = lengths[k]
        vals[i] = unvec(field[k], vector[idx:idx+n-1])
        idx += n
    end
    
    return NamedTuple{names}(Tuple(vals))
end

"""
    unravel(cp, vector) -> (du, dq, m)

Reconstruct boundary condition, source, and mass fraction NamedTuples from a flat control vector.
"""
function unravel(cp::ControlParameters, vector::Vector)
    # Calculate section boundaries
    du_len = total_length(cp.du_lengths)
    dq_len = total_length(cp.dq_lengths)
    m_len = total_length(cp.m_lengths)
    
    idx = 1
    du = nothing; dq = nothing; m = nothing
    
    if !isnothing(cp.du)
        du = unravel_field(cp.du, cp.du_lengths, vector[idx:idx+du_len-1])
        idx += du_len
    end
    
    if !isnothing(cp.dq)
        dq = unravel_field(cp.dq, cp.dq_lengths, vector[idx:idx+dq_len-1])
        idx += dq_len
    end
    
    if !isnothing(cp.m)
        m = unravel_field(cp.m, cp.m_lengths, vector[idx:idx+m_len-1])
    end
    
    return du, dq, m
end
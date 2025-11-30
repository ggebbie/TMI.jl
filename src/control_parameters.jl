"""
    check_shared_references(nt::NamedTuple, name::String)

Check if any non-nothing fields in a NamedTuple share the same memory address.
Throws an error if shared references are detected.

# Detects shared objects so in-place updates on one field do not silently
# mutate another.

# Arguments
- `nt`: NamedTuple to check
- `name`: Name of the variable (for error message)

# Examples
```julia
b = BoundaryCondition(...)
u₀ = (c = b, c_q = b)  # Both reference same object
check_shared_references(u₀, "u₀")  # Throws error
```
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
                error("Shared reference detected in $name: fields $(join(shared_keys, ", ")) " *
                      "point to the same object in memory. Use deepcopy to create independent copies.")
            end
        end
    end
end
check_shared_references(x, name::String) = nothing  # Non-NamedTuple, no check needed

"""
    ControlParameters

Container for optimization control variables (boundary conditions, sources, mass fractions) with vectorization support.
Includes optional first guess values and weighting matrices for inverse problems.

# Fields
- `u`: NamedTuple of boundary condition perturbations (BoundaryCondition objects)
- `q`: NamedTuple of source perturbations (Source objects)
- `m`: NamedTuple of mass fractions (MassFraction objects)
- `u₀`: NamedTuple of first guess boundary conditions
- `q₀`: NamedTuple of first guess sources
- `m₀`: NamedTuple of first guess mass fractions
- `Qᵤ`: NamedTuple of inverse covariance matrices (Symmetric or Diagonal) for boundary condition errors
- `Qₛ`: NamedTuple of inverse covariance matrices (Symmetric or Diagonal) for source errors
- `Qₘ`: Single inverse covariance matrix (Symmetric or Diagonal) for mass fraction errors
- `vector`: Flattened vector representation of all control parameters (if store_vector=true)
- `u_length`: Total length of vectorized `u` (if store_lengths=true)
- `q_length`: Total length of vectorized `q` (if store_lengths=true)
- `m_length`: Total length of vectorized `m` (if store_lengths=true)

# Constructor
```julia
ControlParameters(; u = nothing, q = nothing, m = nothing,
                   u₀ = nothing, q₀ = nothing, m₀ = nothing,
                   Qᵤ = nothing, Qₛ = nothing, Qₘ = nothing,
                   store_vector = true, store_lengths = true)
```
"""
struct ControlParameters
    u::Union{NamedTuple, Nothing}
    q::Union{NamedTuple, Nothing}
    m::Union{NamedTuple, Nothing}
    u₀::Union{NamedTuple, Nothing}
    q₀::Union{NamedTuple, Nothing}
    m₀::Union{NamedTuple, Nothing}
    Qᵤ::Union{NamedTuple, Nothing}
    Qₛ::Union{NamedTuple, Nothing}
    Qₘ::Union{Symmetric, Diagonal, Nothing}
    vector::Union{Vector, Nothing}
    u_length::Union{Int, Nothing}
    q_length::Union{Int, Nothing}
    m_length::Union{Int, Nothing}

    function ControlParameters(; u = nothing, q = nothing, m = nothing,
                                u₀ = nothing, q₀ = nothing, m₀ = nothing,
                                Qᵤ = nothing, Qₛ = nothing, Qₘ = nothing,
                                store_vector = true, store_lengths = true)
        # Ensure control entries and priors are independent objects.
        check_shared_references(u, "u")
        check_shared_references(q, "q")
        check_shared_references(m, "m")
        check_shared_references(u₀, "u₀")
        check_shared_references(q₀, "q₀")
        check_shared_references(m₀, "m₀")
        check_shared_references(Qᵤ, "Qᵤ")
        check_shared_references(Qₛ, "Qₛ")

        # Optionally cache flattened controls for Optim/ForwardDiff style interfaces.
        vector = store_vector ? vectorize_controls(u, q, m) : nothing

        # Segment lengths let `unvec` slice an incoming control vector without recomputing sizes.
        u_length = store_lengths && !isnothing(u) ? length(vec(u)) : nothing
        q_length = store_lengths && !isnothing(q) ? length(vec(q)) : nothing
        m_length = store_lengths && !isnothing(m) ? length(vec(m)) : nothing

        new(u, q, m, u₀, q₀, m₀, Qᵤ, Qₛ, Qₘ, vector, u_length, q_length, m_length)
    end
end

"""
    vectorize_controls(u, q, m) -> Vector

Flatten control parameters into a single vector by concatenating vectorized non-nothing fields.
"""
function vectorize_controls(u, q, m)
    return vcat(vec.([u, q, m])...)
end

"""
    vectorize_controls(cp::ControlParameters) -> Vector

Extract or compute the flattened vector representation from a ControlParameters object.
If the vector is already stored, return it. Otherwise, compute it from the control fields.
"""
function vec(cp::ControlParameters)
    if !isnothing(cp.vector)
        return cp.vector
    else
        return vectorize_controls(cp.u, cp.q, cp.m)
    end
end

"""
    unvec(cp, vector) -> (u, q, m)

Reconstruct boundary condition, source, and mass fraction NamedTuples from a flat control vector.
"""
function unvec(cp::ControlParameters, vector::AbstractVector)
    # Calculate boundaries
    u_len = cp.u_length
    q_len = cp.q_length
    m_len = cp.m_length
    
    idx = 1
    u = nothing; q = nothing; m = nothing
    
    if !isnothing(cp.u)
        u = unvec(cp.u, vector[idx:idx+u_len-1])
        idx += u_len
    end
    
    if !isnothing(cp.q)
        q = unvec(cp.q, vector[idx:idx+q_len-1])
        idx += q_len
    end
    
    if !isnothing(cp.m)
        m = unvec(cp.m, vector[idx:idx+m_len-1])
    end

    return u, q, m
end

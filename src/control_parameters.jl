"""
    check_shared_references(nt::NamedTuple, name::String)

Check if any non-nothing fields in a NamedTuple share the same memory address.
Throws an error if shared references are detected.

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
- `du`: NamedTuple of boundary condition perturbations (BoundaryCondition objects)
- `dq`: NamedTuple of source perturbations (Source objects)
- `m`: NamedTuple of mass fractions (MassFraction objects)
- `u₀`: NamedTuple of first guess boundary conditions
- `q₀`: NamedTuple of first guess sources
- `m₀`: NamedTuple of first guess mass fractions
- `Qⁱᵤ`: NamedTuple of inverse covariance matrices (Symmetric or Diagonal) for boundary condition errors
- `Qⁱₛ`: NamedTuple of inverse covariance matrices (Symmetric or Diagonal) for source errors
- `Qⁱₘ`: Single inverse covariance matrix (Symmetric or Diagonal) for mass fraction errors
- `vector`: Flattened vector representation of all control parameters (if store_vector=true)
- `du_length`: Total length of vectorized `du` (if store_lengths=true)
- `dq_length`: Total length of vectorized `dq` (if store_lengths=true)
- `m_length`: Total length of vectorized `m` (if store_lengths=true)

# Constructor
```julia
ControlParameters(; du = nothing, dq = nothing, m = nothing,
                   u₀ = nothing, q₀ = nothing, m₀ = nothing,
                   Qⁱᵤ = nothing, Qⁱₛ = nothing, Qⁱₘ = nothing,
                   store_vector = true, store_lengths = true)
```

All parameters are optional and default to `nothing`. By default, `store_vector=true` stores
a flattened vector representation of the control parameters. By default, `store_lengths=true`
enables unraveling vectors back to structured NamedTuples using the `unravel` function.
Set either to `false` to reduce memory usage if not needed.

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
struct ControlParameters
    du::Union{NamedTuple, Nothing}
    dq::Union{NamedTuple, Nothing}
    m::Union{NamedTuple, Nothing}
    u₀::Union{NamedTuple, Nothing}
    q₀::Union{NamedTuple, Nothing}
    m₀::Union{NamedTuple, Nothing}
    Qⁱᵤ::Union{NamedTuple, Nothing}
    Qⁱₛ::Union{NamedTuple, Nothing}
    Qⁱₘ::Union{Symmetric, Diagonal, Nothing}
    vector::Union{Vector, Nothing}
    du_length::Union{Int, Nothing}
    dq_length::Union{Int, Nothing}
    m_length::Union{Int, Nothing}

    function ControlParameters(; du = nothing, dq = nothing, m = nothing,
                                u₀ = nothing, q₀ = nothing, m₀ = nothing,
                                Qⁱᵤ = nothing, Qⁱₛ = nothing, Qⁱₘ = nothing,
                                store_vector = true, store_lengths = true)
        # Check for aliasing issues - fields should not share memory
        check_shared_references(du, "du")
        check_shared_references(dq, "dq")
        check_shared_references(m, "m")
        check_shared_references(u₀, "u₀")
        check_shared_references(q₀, "q₀")
        check_shared_references(m₀, "m₀")
        check_shared_references(Qⁱᵤ, "Qⁱᵤ")
        check_shared_references(Qⁱₛ, "Qⁱₛ")

        if store_vector
            vector = vectorize_controls(du, dq, m)
        else
            vector = nothing
        end

        if store_lengths
            du_length = isnothing(du) ? nothing : length(vec(du))
            dq_length = isnothing(dq) ? nothing : length(vec(dq))
            m_length = isnothing(m) ? nothing : length(vec(m))
        else
            du_length = nothing
            dq_length = nothing
            m_length = nothing
        end

        new(du, dq, m, u₀, q₀, m₀, Qⁱᵤ, Qⁱₛ, Qⁱₘ, vector, du_length, dq_length, m_length)
    end
end

"""
    vectorize_controls(du, dq, m) -> Vector

Flatten control parameters into a single vector by concatenating vectorized non-nothing fields.
"""
function vectorize_controls(du, dq, m)
    return vcat(vec.([du, dq, m])...)
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
        return vectorize_controls(cp.du, cp.dq, cp.m)
    end
end

"""
    unvec(cp, vector) -> (du, dq, m)

Reconstruct boundary condition, source, and mass fraction NamedTuples from a flat control vector.
"""
function unvec(cp::ControlParameters, vector::AbstractVector)
    # Calculate boundaries
    du_len = cp.du_length
    dq_len = cp.dq_length
    m_len = cp.m_length
    
    idx = 1
    du = nothing; dq = nothing; m = nothing
    
    if !isnothing(cp.du)
        du = unvec(cp.du, vector[idx:idx+du_len-1])
        idx += du_len
    end
    
    if !isnothing(cp.dq)
        dq = unvec(cp.dq, vector[idx:idx+dq_len-1])
        idx += dq_len
    end
    
    if !isnothing(cp.m)
        m = unvec(cp.m, vector[idx:idx+m_len-1])
    end

    return du, dq, m
end
"""
    ControlParameters{u_names, q_names, m_names}

Container for optimization control variables (boundary conditions, sources, mass fractions) with vectorization support.
"""
struct ControlParameters{u_names, q_names, m_names}
    du::Union{NamedTuple{u_names, <:Tuple{Vararg{BoundaryCondition}}}, Nothing}
    dq::Union{NamedTuple{q_names, <:Tuple{Vararg{Source}}}, Nothing}
    m::Union{NamedTuple{m_names, <:Tuple{Vararg{MassFraction}}}, Nothing}
    vector::Vector
    du_lengths::Union{Dict{Symbol, Int}, Nothing}
    dq_lengths::Union{Dict{Symbol, Int}, Nothing}
    m_lengths::Union{Dict{Symbol, Int}, Nothing}
    
    function ControlParameters(; du = nothing, dq = nothing, m = nothing, store_lengths = false)
        vector = vcat(vec.(filter(!isnothing, [du, dq, m]))...)
        
        u_names = isnothing(du) ? () : typeof(du).parameters[1]
        q_names = isnothing(dq) ? () : typeof(dq).parameters[1]
        m_names = isnothing(m) ? () : typeof(m).parameters[1]
        
        if store_lengths
            du_lengths = isnothing(du) ? 0 : Dict(k => length(vec(v)) for (k,v) in pairs(du))
            dq_lengths = isnothing(dq) ? 0 : Dict(k => length(vec(v)) for (k,v) in pairs(dq))
            m_lengths = isnothing(m) ? 0 : Dict(k => length(vec(v)) for (k,v) in pairs(m))
        else
            du_lengths = nothing
            dq_lengths = nothing
            m_lengths = nothing
        end
        
        new{u_names, q_names, m_names}(du, dq, m, vector, du_lengths, dq_lengths, m_lengths)
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
#A struct that stores and vectorizes control parameters
struct ControlParameters{u_names, q_names, m_names}
    du::Union{NamedTuple{u_names, <:Tuple{Vararg{BoundaryCondition}}}, Nothing}
    dq::Union{NamedTuple{q_names, <:Tuple{Vararg{Source}}}, Nothing}
    m::Union{NamedTuple{m_names, <:Tuple{Vararg{MassFraction}}}, Nothing}
    vector::Vector
    
    function ControlParameters(; du = nothing, dq = nothing, m = nothing)
        vector = vcat(vec.(filter(!isnothing, [du, dq, m]))...)
        
        u_names = isnothing(du) ? () : typeof(u).parameters[1]
        q_names = isnothing(dq) ? () : typeof(q).parameters[1]
        m_names = isnothing(m) ? () : typeof(m).parameters[1]
        
        new{u_names, q_names, m_names}(du, dq, m, vector)
    end
end




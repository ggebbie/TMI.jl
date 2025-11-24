
function prior_boundary_cost(duvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + du) - u₀)' S ((u₀ + u) - u₀)
    return duvec' * Qⁱᵤ * duvec #a shortcut since (u₀ + u) - u₀ = u
end

function gprior_boundary_cost(duvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    gJ = 1 #should always be 1
    gduvec = 2 * Qⁱᵤ * duvec * gJ 
    return gduvec
end

function prior_source_cost(dqvec::Vector,q₀::Union{Source,NamedTuple},Qⁱₛ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return dqvec' * Qⁱₛ * dqvec #a shortcut since (u₀ + u) - u₀ = u
end

function gprior_source_cost(dqvec::Vector,q₀::Union{Source,NamedTuple},Qⁱₛ::Symmetric)
    gJ = 1 #should always be 1
    gqvec = 2 * Qⁱₛ * dqvec * gJ 
    return gqvec
end

function model_observation_cost(n::Union{Vector, Field},Wⁱ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return n ⋅ (Wⁱ * n) # dot product, gives total misfit

end

function gmodel_observation_cost(n::Vector,Wⁱ::Symmetric)
    gJ = 1 #should always be 1
    gnvec = 2 * Wⁱ * n * gJ 
    return gnvec
end

function steadyinversion_residual(c::Field, c_obs::Union{Vector, Field}, 
                                  locs::{Nothing, Vector{Tuple{T,T,T}}}, 
                                  γ::Grid)::Vector{T} where T <: Real
    # n = (c - c_obs)
    if c_obs isa Field
        n = c - c_obs
        return n
    # n = (Ec - c_obs)
    elseif c_obs isa Vector
        #observe recalculates weights everytime its called. could be made more efficient by creating 
        #a "locations" object that caches the interpolation weights.
        #i.e., compute and store weights when locs is initalized
        y = observe(c, locs, γ)
        n = y .- c_obs
        return n
    end
end

function gsteadyinversion_residual(dn::Union{Vector, Field},
                                    c::Field,
                                  locs::Vector{Tuple{T,T,T}}, 
                                  γ::Grid)::Vector{T} where T <: Real
    if dn isa Field
        gc = dn
        return gc
    elseif dn isa Vector
        gy = dn 
        gc = gobserve(gy,c,locs) 
        return gc
    end
end

function global_cost_function(control_vector::Vector{T}, 
                              control_vector_struct::ControlParameters, 
                              u₀::NamedTuple{u₀_names, <:Tuple{Vararg{BoundaryCondition}}}, 
                              q₀::NamedTuple{q₀_names, <:Tuple{Vararg{Source}}}, 
                              m₀::NamedTuple{m₀_names, <:Tuple{Vararg{MassFraction}}},
                              Qⁱᵤ::NamedTuple{u₀_names, <:Tuple{Vararg{Symmetric}}},
                              Qⁱₛ::NamedTuple{q₀_names, <:Tuple{Vararg{Symmetric}}},
                              Wⁱ::NamedTuple{cobs_names, <:Tuple{Vararg{Symmetric}}},
                              cobs::NamedTuple{cobs_names, <:Tuple{Vararg{Union{Vector, Matrix}}}},
                              γ::Grid, locs = nothing) where {T <: Real, u₀_names, q₀_names, m₀_names, cobs_names}
    m, du, dq = unravel(control_vector_struct, control_vector)
    mvec = vec(m)
    J, dJdcontrols = global_cost_function(mvec, du, dq, 
                                                 u₀, q₀, m₀,
                                                 Qⁱᵤ, Qⁱₛ, Wⁱ,
                                                 cobs, γ, locs)

    return J, dJdcontrols.vector
end


function global_cost_function(mvec::Vector{T}, 
                              du::NamedTuple{du_names, <:Tuple{Vararg{BoundaryCondition}}}, 
                              dq::NamedTuple{dq_names, <:Tuple{Vararg{Source}}}, 
                              u₀::NamedTuple{u₀_names, <:Tuple{Vararg{BoundaryCondition}}}, 
                              q₀::NamedTuple{q₀_names, <:Tuple{Vararg{Source}}}, 
                              m₀::NamedTuple{m₀_names, <:Tuple{Vararg{MassFraction}}},
                              Qⁱᵤ::NamedTuple{u₀_names, <:Tuple{Vararg{Symmetric}}},
                              Qⁱₛ::NamedTuple{q₀_names, <:Tuple{Vararg{Symmetric}}},
                              Wⁱ::NamedTuple{cobs_names, <:Tuple{Vararg{Symmetric}}},
                              cobs::NamedTuple{cobs_names, <:Tuple{Vararg{Union{Vector, Matrix}}}},
                              γ::Grid, locs = nothing) where {T <: Real, du_names, u₀_names, dq_names, q₀_names, m₀_names, cobs_names}

    m = unvec(m₀, mvec)
    duvec = vec(du)
    dqvec = vec(dq)

    b = deepcopy(u₀)
    q = deepcopy(q₀)
    for key in du_names #only update if in the control vector 
        adjustboundarycondition!(b[key],du[key]) #b += u # easy case where u and b are on the same boundary
    end

    for key in dq_names
        adjustsource!(q[key],dq[key]) #b += u # easy case where u and b are on the same boundary
    end

    A = watermassmatrix(m)
    Alu = lu(A)

    #need to check this works for named tuples
    c = steadyinversion(Alu,b,γ,q=q,r=r) 
    n = steadyinversion_residual(c, cobs, locs, γ)

    Jn = model_observation_cost(n,Wⁱ)
    Ju = prior_boundary_cost(duvec,u₀, Qⁱᵤ)
    Jq = prior_boundary_cost(dqvec,u₀, Qⁱₛ)

    #cost function need to add the lagrangrian mulitpliers as well 
    # i think this is partially given by tracer contribution 
    J = Jn + Ju + Jq

    #need to loop and accumulate here
    gduvec = gprior_boundary_cost(duvec,u₀,Qⁱᵤ)
    gdqvec = gprior_source_cost(dqvec,u₀,Qⁱₛ)

    gn = gmodel_observation_cost(n,Wⁱ)
    gc = gsteadyinversion_residual(gn, c, locs, γ)
    #needs to accumulate (I think)
    gdu_bc,gdq_src, gA = gsteadyinversion(gc,c,A, Alu, b,γ;q=q,r=1.0)

    #a quick way to subset the correct controls
    gduvec .+= vec(NamedTuple{Tuple(du_names)}(gdu_bc))
    gdqvec .+= vec(NamedTuple{Tuple(du_names)}(gdq_src))

    gm = gwatermassmatrix(gA, m, γ)

    dJdcontrols = ControlParameters(; du = unvec(du, gduvec), dq = unvec(dq, gdqvec), m = gm)

    return J, dJdcontrols
end
"""
    prior_mass_fraction_cost(mvec, m₀, Qⁱₘ) -> Real

Compute the quadratic cost of mass fraction deviations from the prior: `(m - m₀)' Qⁱₘ (m - m₀)`.
"""
function prior_mass_fraction_cost(mvec::Vector,m₀::NamedTuple,Qⁱₘ::Symmetric)
    m₀_vec = vec(m₀)
    Δm = (mvec .- m₀_vec)
    return Δm' * Qⁱₘ * Δm
end

"""
    gprior_mass_fraction_cost(mvec, m₀, Qⁱₘ) -> Vector

Compute the gradient of the mass fraction prior cost: `∂/∂m = 2 Qⁱₘ (m - m₀)`.
"""
function gprior_mass_fraction_cost(mvec::Vector,m₀::NamedTuple,Qⁱₘ::Symmetric)
    m₀_vec = vec(m₀)
    Δm = (mvec .- m₀_vec)
    gmvec = 2 * Qⁱₘ * Δm
    return gmvec
end


"""
    prior_boundary_cost(duvec, u₀, Qⁱᵤ) -> Real

Compute the quadratic cost of boundary condition perturbations: `du' Qⁱᵤ du`.
"""
function prior_boundary_cost(duvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    return duvec' * Qⁱᵤ * duvec
end

"""
    gprior_boundary_cost(duvec, u₀, Qⁱᵤ) -> Vector

Compute the gradient of the boundary condition prior cost: `∂/∂du = 2 Qⁱᵤ du`.
"""
function gprior_boundary_cost(duvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    gduvec = 2 * Qⁱᵤ * duvec
    return gduvec
end

"""
    prior_source_cost(dqvec, q₀, Qⁱₛ) -> Real

Compute the quadratic cost of source perturbations: `dq' Qⁱₛ dq`.
"""
function prior_source_cost(dqvec::Vector,q₀::Union{Source,NamedTuple},Qⁱₛ::Symmetric)
    return dqvec' * Qⁱₛ * dqvec
end

"""
    gprior_source_cost(dqvec, q₀, Qⁱₛ) -> Vector

Compute the gradient of the source prior cost: `∂/∂dq = 2 Qⁱₛ dq`.
"""
function gprior_source_cost(dqvec::Vector,q₀::Union{Source,NamedTuple},Qⁱₛ::Symmetric)
    gqvec = 2 * Qⁱₛ * dqvec
    return gqvec
end

"""
    model_observation_cost(n, Wⁱ) -> Real

Compute the weighted observation misfit cost: `n' Wⁱ n`.
"""
function model_observation_cost(n::Union{Vector, Field},Wⁱ::Symmetric)
    return n ⋅ (Wⁱ * n)
end

"""
    gmodel_observation_cost(n, Wⁱ) -> Vector

Compute the gradient of the observation cost: `∂/∂n = 2 Wⁱ n`.
"""
function gmodel_observation_cost(n::Vector,Wⁱ::Symmetric)
    gnvec = 2 * Wⁱ * n
    return gnvec
end

"""
    steadyinversion_residual(c, c_obs, locs, γ) -> Vector

Compute the residual between modeled tracer field `c` and observations `c_obs`.
"""
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
    # elseif c_obs isa Observations
        # y = observe(c, c_obs.wis, γ)
        # n = y .- c_obs.values
        # return n
    end
end

"""
    gsteadyinversion_residual(dn, c, locs, γ) -> Vector

Compute the gradient of the residual with respect to the tracer field `c`.
"""
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

"""
    global_cost_function(control_vector, control_vector_struct, u₀, q₀, m₀, Qⁱᵤ, Qⁱₛ, Wⁱ, cobs, γ, locs) -> (J, gradient)

Evaluate total cost and gradient from a flat control vector; wrapper that unravels controls first.
"""
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
    du, dq, m = unravel(control_vector_struct, control_vector)
    mvec = vec(m)
    J, dJdcontrols = global_cost_function(mvec, du, dq, 
                                                 u₀, q₀, m₀,
                                                 Qⁱᵤ, Qⁱₛ, Wⁱ,
                                                 cobs, γ, locs)

    return J, dJdcontrols.vector
end


"""
    global_cost_function(mvec, du, dq, u₀, q₀, m₀, Qⁱᵤ, Qⁱₛ, Wⁱ, cobs, γ, locs) -> (J, ControlParameters)

Compute the total cost J (observation + prior terms) and its gradient with respect to all controls.
"""
function global_cost_function(mvec::Vector{T},
                              du::NamedTuple{du_names, <:Tuple{Vararg{BoundaryCondition}}},
                              dq::NamedTuple{dq_names, <:Tuple{Vararg{Source}}},
                              u₀::NamedTuple{u₀_names, <:Tuple{Vararg{BoundaryCondition}}},
                              q₀::NamedTuple{q₀_names, <:Tuple{Vararg{Source}}},
                              m₀::NamedTuple{m₀_names, <:Tuple{Vararg{MassFraction}}},
                              Qⁱᵤ::NamedTuple{u₀_names, <:Tuple{Vararg{Symmetric}}},
                              Qⁱₛ::NamedTuple{q₀_names, <:Tuple{Vararg{Symmetric}}},
                              Wⁱ::NamedTuple{cobs_names, <:Tuple{Vararg{Symmetric}}},
                              cobs::NamedTuple{cobs_names, <:Tuple{Vararg{Union{Vector, Field}}}},
                              γ::Grid, locs = nothing) where {T <: Real, du_names, u₀_names, dq_names, q₀_names, m₀_names, cobs_names}

    m = unvec(m₀, mvec)
    duvec = vec(du) #vectorize boundary condition adjustments u = u_0 + du
    dqvec = vec(dq) #vectorize source adjustments q = q_0 + dq

    b = deepcopy(u₀)
    q = deepcopy(q₀)

    #might already exist! 
    for key in du_names #only update if in the control vector 
        adjustboundarycondition!(b[key],du[key]) #b += u # easy case where u and b are on the same boundary
    end

    for key in dq_names
        adjustsource!(q[key],dq[key]) #b += u # easy case where u and b are on the same boundary
    end

    A = watermassmatrix(m)
    Alu = lu(A)

    #need to check this works for named tuples
    #it should now. 
    c = steadyinversion(Alu,b,γ,q=q,r=r) 
    n = steadyinversion_residual(c, cobs, locs, γ)

    Jn = model_observation_cost(n,Wⁱ)
    Ju = prior_boundary_cost(duvec,u₀, Qⁱᵤ)
    Jq = prior_source_cost(dqvec,q₀, Qⁱₛ)
    # Jm = prior_mass_fraction_cost(mvec,m₀,Qⁱₘ)
    #Jλ = λ^Τ * F . I think F is tracer_contributions - q. need to double checl

    #TO DO: add mass fraction prior cost, lagrange multiplier cost 
    #and tracer constraints (i.e., freezing point of seawater, gravitational stability)
    J = Jn + Ju + Jq #+ (λ^Τ * F)

    #need to loop and accumulate here
    gduvec = gprior_boundary_cost(duvec,u₀,Qⁱᵤ)
    gdqvec = gprior_source_cost(dqvec,q₀,Qⁱₛ)

    gn = gmodel_observation_cost(n,Wⁱ)
    gc = gsteadyinversion_residual(gn, c, locs, γ)
    #needs to accumulate (I think)
    gdu_bc,gdq_src, gA = gsteadyinversion(gc,c,A, Alu, b,γ;q=q,r=1.0)

    #a quick way to subset the correct controls
    gduvec .+= vec(NamedTuple{Tuple(du_names)}(gdu_bc))
    gdqvec .+= vec(NamedTuple{Tuple(du_names)}(gdq_src))

    gm = gwatermassmatrix(gA, m, γ)

    # 
    
    dJdcontrols = ControlParameters(; du = unvec(du, gduvec), dq = unvec(dq, gdqvec), m = gm)

    return J, dJdcontrols
end
"""
    prior_mass_fraction_cost(mvec, m₀, Qⁱₘ) -> Real

Compute the quadratic cost of mass fraction deviations from the prior: `(m - m₀)' Qⁱₘ (m - m₀)`.
"""
function prior_mass_fraction_cost(m::Vector,m₀::Vector,Qⁱₘ::Union{Diagonal, Symmetric})
    Δm = (m .- m₀)
    return Δm' * Qⁱₘ * Δm
end

function prior_mass_fraction_cost(m::NamedTuple,m₀::NamedTuple,Qⁱₘ::Union{Diagonal, Symmetric})
    return prior_mass_fraction_cost(vec(m), vec(m₀), Qⁱₘ)
end


"""
    gprior_mass_fraction_cost(mvec, m₀, Qⁱₘ) -> Vector

Compute the gradient of the mass fraction prior cost: `∂/∂m = 2 Qⁱₘ (m - m₀)`.
"""
function gprior_mass_fraction_cost(m::Vector,m₀::Vector,Qₘ::Union{Diagonal, Symmetric})
    Δm = (m .- m₀)
    gmvec = 2 * Qₘ * Δm
    return gmvec
end

function gprior_mass_fraction_cost(m::NamedTuple,m₀::NamedTuple,Qₘ::Union{Diagonal, Symmetric})
    gmvec = gprior_mass_fraction_cost(vec(m), vec(m₀), Qₘ)
    return unvec(m₀, gmvec)
end


"""
    prior_boundary_cost(duvec, u₀, Qⁱᵤ) -> Real

Compute the quadratic cost of boundary condition perturbations: `du' Qⁱᵤ du`.
"""
function prior_boundary_cost(duvec::AbstractVector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Union{Diagonal, Symmetric})
    return duvec' * Qⁱᵤ * duvec
end

function prior_boundary_cost(
      duvec::Vector,
      u₀::NamedTuple,
      Qᵤ::NamedTuple{tracer_names,T}
  ) where {tracer_names, T}
    J = 0
    lo = 1
    for name in tracer_names
        len = size(Qᵤ[name], 1)
        hi = lo + len - 1
        du_slice = @view duvec[lo:hi]
        J +=  prior_boundary_cost(du_slice, u₀[name], Qᵤ[name])
        lo = hi + 1
    end
    return J
end


# function prior_boundary_cost(du::NamedTuple,u₀::NamedTuple,Qᵤ::NamedTuple)
#     tracer_names = keys(du)
#     J = 0
#     for (i, name) in enumerate(tracer_names)
#         J += prior_boundary_cost(du[name], u₀[name], Qᵤ[name])
#     end
#     return J
# end



"""
    gprior_boundary_cost(duvec, u₀, Qⁱᵤ) -> Vector

Compute the gradient of the boundary condition prior cost: `∂/∂du = 2 Qⁱᵤ du`.
"""
function gprior_boundary_cost(duvec::AbstractVector,
                              u₀::Union{BoundaryCondition,NamedTuple},
                              Qⁱᵤ::Union{Diagonal, Symmetric})
    return 2 * Qⁱᵤ * duvec
end

function gprior_boundary_cost(
      duvec::AbstractVector,
      u₀::NamedTuple{names},
      Qᵤ::NamedTuple{names,<:Tuple{Vararg{Union{Symmetric,Diagonal}}}}
  ) where {names}
    gduvec = similar(duvec)
    lo = 1
    for name in names
        len = size(Qᵤ[name], 1)
        hi = lo + len - 1
        gdu_slice  = @view gduvec[lo:hi] #using @view for a performance boost
        du_slice = @view duvec[lo:hi]
        gdu_slice .= gprior_boundary_cost(du_slice, u₀[name], Qᵤ[name])
        lo = hi + 1
    end
    return gduvec
end

function gprior_boundary_cost(du::NamedTuple, u₀::NamedTuple, Qᵤ::NamedTuple)
    tracer_names = keys(du)
    n_tracers = length(tracer_names)
    gdu_results = Vector(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        gdu_i = gprior_boundary_cost(vec(du[name]), u₀[name], Qᵤ[name])
        gdu_results[i] = unvec(u₀[name], gdu_i)
    end
    gdu_nt = (; zip(tracer_names, gdu_results)...)
    return gdu_nt

end


"""
    prior_source_cost(dqvec, q₀, Qⁱₛ) -> Real

Compute the quadratic cost of source perturbations: `dq' Qⁱₛ dq`.
"""
function prior_source_cost(dqvec::Union{Vector, SubArray},q₀::Union{Source,NamedTuple},Qⁱₛ::Union{Diagonal, Symmetric})
    return dqvec' * Qⁱₛ * dqvec
end

function prior_source_cost(
      dqvec::Vector,
      q₀::NamedTuple,
      Qₛ::NamedTuple{tracer_names,<:Tuple{Vararg{Union{Nothing, Symmetric,Diagonal}}}}
  ) where {tracer_names}
    J = 0
    lo = 1
    for name in tracer_names
        if !isnothing(q₀[name])
            len = size(Qₛ[name], 1)
            hi = lo + len - 1
            dq_slice = @view dqvec[lo:hi]
            J +=  prior_source_cost(dq_slice, q₀[name], Qₛ[name])
            lo = hi + 1
        end
    end
    return J
end

function prior_boundary_cost(dq::NamedTuple,q₀::NamedTuple,Qₛ::NamedTuple)
    tracer_names = keys(dq)
    J = 0
    for (i, name) in enumerate(tracer_names)
        if !isnothing(q₀[name])
            J += prior_source_cost(dq[name], q₀[name], Qₛ[name])
        end
    end
    return J
end

"""
    gprior_source_cost(dqvec, q₀, Qⁱₛ) -> Vector

Compute the gradient of the source prior cost: `∂/∂dq = 2 Qⁱₛ dq`.
"""
function gprior_source_cost(dqvec::Union{Vector, SubArray},q₀::Union{Source,NamedTuple},Qⁱₛ::Union{Diagonal, Symmetric})
    gdqvec = 2 * Qⁱₛ * dqvec
    return gdqvec
end

function gprior_source_cost(
      dqvec::Vector,
      q₀::NamedTuple,
      Qₛ::NamedTuple{tracer_names,T}
  ) where {tracer_names, T}
    gdqvec = similar(dqvec)
    lo = 1
    for name in tracer_names
        if !isnothing(Qₛ[name])
            len = size(Qₛ[name], 1)
            hi = lo + len - 1
            gdq_slice  = @view gdqvec[lo:hi] #using @view for a performance boost
            dq_slice = @view dqvec[lo:hi]
            gdq_slice .= gprior_source_cost(dq_slice, q₀[name], Qₛ[name])
            lo = hi + 1
        end
    end
    return gdqvec
end

function gprior_source_cost(
    dq::NamedTuple, 
    q₀::NamedTuple, 
    Qₛ::NamedTuple{tracer_names,T}
    ) where {tracer_names, T}
    n_tracers = length(tracer_names)
    gdq_results = Vector(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        if !isnothing(Qₛ[name])
            gdq_i = gprior_source_cost(vec(dq[name]), q₀[name], Qₛ[name])
            gdq_results[i] = unvec(q₀[name], gdq_i)
        else
            gdq_results[i] = nothing
        end
    end
    gdq_nt = (; zip(tracer_names, gdq_results)...)
    return gdq_nt
end

"""
    model_observation_cost(n, Wⁱ) -> Real

Compute the weighted observation misfit cost: `n' Wⁱ n`.
"""
function model_observation_cost(n::Union{Vector, Field},Wⁱ::Union{Diagonal, Symmetric})
    return n ⋅ (Wⁱ * n)
end

function model_observation_cost(n::NamedTuple,c_obs::NamedTuple)
    tracer_names = keys(n)
    J = 0
    for (i, name) in enumerate(tracer_names)
        J += model_observation_cost(n[name], c_obs[name].W)
    end
    return J
end

"""
    gmodel_observation_cost(n, Wⁱ) -> Vector

Compute the gradient of the observation cost: `∂/∂n = 2 Wⁱ n`.
"""
function gmodel_observation_cost(n::Union{Vector, Field},Wⁱ::Union{Diagonal, Symmetric})
    gnvec = 2 * Wⁱ * n
    return gnvec
end

function gmodel_observation_cost(n::NamedTuple, c_obs::NamedTuple)
    tracer_names = keys(n)
    n_tracers = length(tracer_names)
    gn_results = Vector(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        gn_results[i] = gmodel_observation_cost(n[name], c_obs[name].W)
    end
    dn_nt = (; zip(tracer_names, gn_results)...)
    return dn_nt
end


"""
    steadyinversion_residual(c, c_obs, locs, γ) -> Vector

Compute the residual between modeled tracer field `c` and observations `c_obs`.
"""
function model_data_misfit(c::Field, c_obs::Union{Vector, Field, Observations},
                                  γ::Grid; locs::Union{Nothing, Vector{G}}=nothing) where G
    # n = (c - c_obs)
    if c_obs isa Field
        n = c - c_obs
        # n = (Ec - c_obs)
    elseif c_obs isa Vector
        y = observe(c, locs, γ)
        n = y .- c_obs
    elseif c_obs isa Observations
        c_obs_val = c_obs.values
        if c_obs_val isa Vector
            y = observe(c, c_obs.wis, γ) #use pre-computed weights
            n = y .- c_obs_val
        else
            n = c - c_obs_val
        end
    end
    return n
end

function model_data_misfit(c::NamedTuple, c_obs::NamedTuple,
                                  γ::Grid; locs::Union{Nothing, NamedTuple}=nothing)

    tracer_names = keys(c_obs)
    n_tracers = length(tracer_names)
    n_results = Vector(undef, n_tracers) #eltype of c_obs or c_obs.values or something

    for (i, name) in enumerate(tracer_names)
        c_i = c[name]
        c_obs_i = c_obs[name]
        n_results[i] = model_data_misfit(c_i, c_obs_i, γ; locs = locs)
    end
    n_nt = NamedTuple{tracer_names}(Tuple(n_results))
    return n_nt

end


"""
    gmodel_data_misfit(dn, c, locs, γ) -> Vector

Compute the gradient of the residual with respect to the tracer field `c`.
"""
function gmodel_data_misfit(dn::Union{Vector, Field},
                                    c::Field, c_obs::Union{Vector, Field, Observations}, 
                                  γ::Grid; locs::Union{Nothing, NamedTuple}=nothing)
    if c_obs isa Field
        gc = dn
    elseif c_obs isa Vector
        gy = dn 
        gc = gobserve(gy,c,locs) 
    elseif c_obs isa Observations
        c_obs_val = c_obs.values
        if c_obs_val isa Field
            gc = dn
        elseif c_obs_val isa Vector
            gy = dn 
            gc = gobserve(gy,c,c_obs.locs; wis = c_obs.wis) 
        end
    end

    return gc

end

function gmodel_data_misfit(gn::NamedTuple, c::NamedTuple, c_obs::NamedTuple,
                                  γ::Grid; locs::Union{Nothing, NamedTuple}=nothing)

    tracer_names = keys(gn)
    n_tracers = length(tracer_names)
    gc_results = Vector(undef, n_tracers)

    for (i, name) in enumerate(tracer_names)
        dn_i = gn[name]
        c_i = c[name]
        c_obs_i = c_obs[name]
        gc_results[i] = gmodel_data_misfit(dn_i, c_i, c_obs_i, γ; locs = locs)
    end
    gc_nt = (; zip(tracer_names, gc_results)...)
    return gc_nt

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

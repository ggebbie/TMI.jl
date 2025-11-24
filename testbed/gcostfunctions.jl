import Pkg; Pkg.activate(".")

using Revise
using TMI
using SparseArrays

import TMI:MassFraction
import TMI:step_cartesian

function prior_boundary_cost(uvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return uvec' * Qⁱᵤ * uvec #a shortcut since (u₀ + u) - u₀ = u
end

function gprior_boundary_cost(uvec::Vector,u₀::Union{BoundaryCondition,NamedTuple},Qⁱᵤ::Symmetric)
    gJ = 1 #should always be 1
    guvec = 2 * Qⁱᵤ * uvec * gJ 
    return guvec
end

function prior_source_cost(qvec::Vector,q₀::Union{BoundaryCondition,NamedTuple},Qⁱₛ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return qvec' * Qⁱₛ * qvec #a shortcut since (u₀ + u) - u₀ = u
end

function gprior_source_cost(qvec::Vector,q₀::Union{BoundaryCondition,NamedTuple},Qⁱₛ::Symmetric)
    gJ = 1 #should always be 1
    gqvec = 2 * Qⁱₛ * qvec * gJ 
    return gqvec
end

function steadyinversion_residual(c::Field, c_obs::Union{Vector, Field}, 
                                  locs::{Nothing, Vector{Tuple{T,T,T}}}, 
                                  γ::Grid)::Vector{T} where T <: Real
    if c_obs isa Field
        r = c - c_obs
        return r
    elseif c_obs isa Vector
        #this recalculates weights everytime its called. could be made more efficient by creating 
        #a "locations" object that caches the interpolation weights 
        y = observe(c, locs, γ)
        r = y .- c_obs
        return r
    end
end


function gsteadyinversion_residual(dr::Union{Vector, Field},
                                    c::Field,
                                  locs::Vector{Tuple{T,T,T}}, 
                                  γ::Grid)::Vector{T} where T <: Real
    if dr isa Field
        gc = dr
        return gc
    elseif dr isa Vector
        gy = dr 
        gc = gobserve(gy,c,locs) 
        return gc
    end
end

function model_observation_cost(r::Union{Vector, Field},Wⁱ::Symmetric)
    #calculate the cost with respect to a prior boundary condition 
    # ((u₀ + u) - u₀)' S ((u₀ + u) - u₀)
    return r ⋅ (Wⁱ * r) # dot product, gives total misfit

end

function gmodel_observation_cost(r::Vector,Wⁱ::Symmetric)
    gJ = 1 #should always be 1
    grvec = 2 * Wⁱ * r * gJ 
    return grvec
end
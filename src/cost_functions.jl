import TMI: Source

"""
    prior_mass_fraction_cost(m, m₀, Qⁱₘ) -> Real

Quadratic penalty for mass-fraction departures from the prior: `(m - m₀)' Qⁱₘ (m - m₀)`.
"""
function prior_mass_fraction_cost(m::Vector,m₀::Vector,Qⁱₘ::Union{Diagonal, Symmetric})
    Δm = (m .- m₀)
    return Δm' * Qⁱₘ * Δm
end

function prior_mass_fraction_cost(m::NamedTuple,m₀::NamedTuple,Qⁱₘ::Union{Diagonal, Symmetric})
    return prior_mass_fraction_cost(vec(m), vec(m₀), Qⁱₘ)
end


"""
    gprior_mass_fraction_cost(m, m₀, Qₘ) -> Vector

Gradient of the mass-fraction prior: `∂/∂m = 2 Qₘ (m - m₀)`.
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
    prior_boundary_cost(du, u₀, Qⁱᵤ) -> Real

Quadratic penalty on boundary-condition perturbations: `du' Qⁱᵤ du`.
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


function prior_boundary_cost(du::NamedTuple,u₀::NamedTuple,Qᵤ::NamedTuple)
    tracer_names = keys(du)
    J = 0
    for (i, name) in enumerate(tracer_names)
        J += prior_boundary_cost(vec(du[name]), u₀[name], Qᵤ[name])
    end
    return J
end



"""
    gprior_boundary_cost(du, u₀, Qⁱᵤ) -> Vector

Gradient of the boundary prior: `∂/∂du = 2 Qⁱᵤ du`.
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
    prior_source_cost(dq, q₀, Qⁱₛ) -> Real

Quadratic penalty on source perturbations: `dq' Qⁱₛ dq`.
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

function prior_source_cost(dq::NamedTuple,q₀::NamedTuple,Qₛ::NamedTuple)
    tracer_names = keys(dq)
    J = 0
    for (i, name) in enumerate(tracer_names)
        if !isnothing(q₀[name])
            J += prior_source_cost(vec(dq[name]), q₀[name], Qₛ[name])
        end
    end
    return J
end

"""
    gprior_source_cost(dq, q₀, Qⁱₛ) -> Vector

Gradient of the source prior: `∂/∂dq = 2 Qⁱₛ dq`.
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

Weighted observation misfit: `n' Wⁱ n`.
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

Gradient of the observation misfit: `∂/∂n = 2 Wⁱ n`.
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
    model_data_misfit(c, c_obs, γ; locs=nothing) -> residual

Residual between model state and observations (fields or sampled vectors); handles Field,
Vector, or Observations inputs.
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
    gmodel_data_misfit(dn, c, c_obs, γ; locs=nothing) -> gradient

Adjoint of `model_data_misfit`: propagate residual sensitivities `dn` back to the state.
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
    unconstrained_global_costfunction(du, dq, m, controls, c_obs, γ; locs=nothing, return_gradients=false)

Objective for the unconstrained inversion: update `b`/`q` with controls, assemble `A`,
solve the steady state, add observation misfit and priors, and optionally return adjoint
gradients w.r.t. `du`, `dq`, and `m` when `return_gradients=true`.
"""
function unconstrained_global_costfunction(control_vector::Vector, 
                                           controls::ControlParameters, c_obs, γ; 
                                           locs = nothing, return_gradients = false)
    u, q, m = unvec(controls, control_vector)
    if !return_gradients
        J, gu, gq, gm = unconstrained_global_costfunction(u, q, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)
        return J
    else
        J, gu, gq, gm = unconstrained_global_costfunction(u, q, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)
        gcontrols = vcat(vec.([gu, gq, gm])...)

        return J, gcontrols
    end
    
end

@inline function _write_gradient!(G, gu, gq, gm, controls::ControlParameters)
    idx = 1
    if !isnothing(controls.u)
        guv = vec(gu)
        len = length(guv)
        copyto!(G, idx, guv, 1, len)
        idx += len
    end
    if !isnothing(controls.q)
        gqv = vec(gq)
        len = length(gqv)
        copyto!(G, idx, gqv, 1, len)
        idx += len
    end
    if !isnothing(controls.m)
        gmv = vec(gm)
        len = length(gmv)
        copyto!(G, idx, gmv, 1, len)
    end
    return nothing
end

function constrained_global_costfunction(control_vector::Vector, 
                                           controls::ControlParameters, c_obs, γ; 
                                           locs = nothing, return_gradients = false)
    u, q, x = unvec(controls, control_vector) #x is a real number vector 
    
    α = 5. #an optional parameter to shrink gradients related to this transformation
    m = softmax_massfractions(x; α = α) #transform x ∈ R to mass fraction that are non-negative and sum to 1

    if !return_gradients
        J, gu, gq, gm = unconstrained_global_costfunction(u, q, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)
        
        return J
    else
        J, gu, gq, gm = unconstrained_global_costfunction(u, q, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)
        gx = gsoftmax_massfractions(gm, m; α = α)

        gcontrols = vcat(vec.([gu, gq, gx])...)

        return J, gcontrols
    end
    
end

"""
    unconstrained_global_costfunction(du, dq, m, controls, c_obs, γ; locs=nothing, return_gradients=false)

Objective for the unconstrained inversion: update `b`/`q` with controls, assemble `A`,
solve the steady state, add observation misfit and priors, and optionally return adjoint
gradients w.r.t. `du`, `dq`, and `m` when `return_gradients=true`.
"""
function unconstrained_global_costfunction(ub, uq, m, controls::ControlParameters, 
                                           c_obs, γ; locs = nothing, return_gradients = false)
    # Forward pass: compute state at current control vector
    b = deepcopy(controls.u₀) #deep-copying to 
    q = deepcopy(controls.q₀)
    du = deepcopy(ub)
    dq = deepcopy(uq)

    # Apply control updates: du,dq store perturbations relative to priors.
    du_names = keys(ub)
    dq_names = keys(uq)

    for key in keys(b)
        if key ∈ du_names
            adjustboundarycondition!(du[key], controls.u₀[key]; r = -1.0) #du = ub - u₀
            adjustboundarycondition!(b[key], du[key]; r = 1.0) #du = ub - u₀
        end
    end

    for key in keys(q)
        if !isnothing(q[key])
            if key ∈ du_names
                adjustsource!(dq[key], controls.q₀[key]; r = -1.0) #dq = uq - q₀
                adjustsource!(q[key], dq[key]; r = 1.0)  #uq = uq + q₀
            end
        end
    end

    A = watermassmatrix(m, γ)
    Alu = lu(A)

    # Solve for each tracer
    c = steadyinversion(Alu,b, q, γ)
    n = model_data_misfit(c, c_obs, γ; locs=locs)
    J = model_observation_cost(n,c_obs) + 
        prior_source_cost(dq, controls.q₀, controls.Qₛ) + 
        prior_boundary_cost(du, controls.u₀, controls.Qᵤ) + 
        prior_mass_fraction_cost(m,controls.m₀,controls.Qₘ)
        
    if !return_gradients 
        return J, nothing, nothing, nothing #just return cost
    else
        # Propagate gradients through adjoint operations
        gdu = zero(du)
        gdq = zero(dq)
        gm = similar(controls.m₀)

        # Prior contributions (no model dynamics).
        gdu_1 = gprior_boundary_cost(du,controls.u₀, controls.Qᵤ)
        gdq_1 = gprior_source_cost(dq,controls.q₀, controls.Qₛ)
        gm_1 = gprior_mass_fraction_cost(m, controls.m₀, controls.Qₘ)

        # Start adjoint by differentiating the observation-cost wrt residuals.
        gn = gmodel_observation_cost(n, c_obs)
        gc = gmodel_data_misfit(gn, c, c_obs, γ; locs=locs)
        gdu_2, gdq_2, gA_total = gsteadyinversion(gc, c, A, Alu, b, q, γ)
        gm_2 = gwatermassmatrix(gA_total, m, γ)

        # Accumulate adjoint contributions into control gradients.
        for key in du_names
            gadjustboundarycondition!(gdu[key], gdu_1[key])
            gadjustboundarycondition!(gdu[key], gdu_2[key])
        end

        for key in dq_names
            if !isnothing(q[key])
                gadjustsource!(gdq[key], gdq_1[key], controls.q₀[key])
                gadjustsource!(gdq[key], gdq_2[key], controls.q₀[key])
            end
        end

        for key in keys(gm)
            gm[key].fraction .= gm_1[key].fraction .+ gm_2[key].fraction
        end

        #since du = u - u₀, ∂J/∂du = ∂J/∂u. Same for dq if q is not on logscale 
        return J, gdu, gdq, gm
    end
end


"""
    optim_fg_unconstrained_global_costfunction!(F, G, control_vector, controls, c_obs, γ; locs=nothing)

Combined objective/gradient for `unconstrained_global_costfunction` for use in Optim.jl.
Call with explicit `controls`, `c_obs`, `γ` (and optional `locs`); computes cost and, when
`G` is provided, the gradient in one pass.
Usage example:
`fg! = (F, G, x) -> optim_fg_unconstrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=locs)`
and pass `Optim.only_fg!(fg!)` to `Optim.optimize`.
"""
function optim_fg_unconstrained_global_costfunction!(F, G, control_vector,
    controls::ControlParameters, c_obs, γ; locs = nothing)
    # Compute objective (and gradient if requested) in one pass.
    u, q, m = unvec(controls, control_vector)
    return_gradients = (G !== nothing) ? true : false

    J, gu, gq, gm = unconstrained_global_costfunction(u, q, m, controls, c_obs, γ; 
                                            locs = locs, return_gradients = return_gradients)

    if G !== nothing
        _write_gradient!(G, gu, gq, gm, controls)
    end

    if F !== nothing
        return J
    end
end


"""
    optim_fg_unconstrained_global_costfunction!(F, G, control_vector, controls, c_obs, γ; locs=nothing)

Combined objective/gradient for `unconstrained_global_costfunction` for use in Optim.jl.
Call with explicit `controls`, `c_obs`, `γ` (and optional `locs`); computes cost and, when
`G` is provided, the gradient in one pass.
Usage example:
`fg! = (F, G, x) -> optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=locs)`
and pass `Optim.only_fg!(fg!)` to `Optim.optimize`.
"""
function optim_fg_constrained_global_costfunction!(F, G, control_vector,
    controls::ControlParameters, c_obs, γ; locs = nothing)
    # Compute objective (and gradient if requested) in one pass.
    u, q, x = unvec(controls, control_vector) #x is a real number vector 
    return_gradients = (G !== nothing) ? true : false
    α = 5. #an optional parameter to shrink gradients related to this transformation
    m = softmax_massfractions(x; α = α) #transform x ∈ R to mass fraction that are non-negative and sum to 1

    J, gu, gq, gm = unconstrained_global_costfunction(u, q, m, controls, c_obs, γ; 
                                              locs = locs, return_gradients = return_gradients)

    if G !== nothing
        gx = gsoftmax_massfractions(gm, m; α = α)
        _write_gradient!(G, gu, gq, gx, controls)
    end

    if F !== nothing
        return J
    end
end

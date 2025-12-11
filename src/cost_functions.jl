
"""
    _write_gradient!(G, gu, gq, gm, controls)

A helper function to flatten the structured gradients (`gu`, `gq`, `gm`) into a single vector `G`. This is used to prepare the gradient for consumption by an optimization algorithm. The function writes the components of the gradient in-place into the vector `G`.
"""
@inline function _write_gradient!(G, gu, gq, gm, controls::Controls)
    idx = 1
    if !isnothing(controls.boundary.ub)
        guv = vec(gu)
        len = length(guv)
        copyto!(G, idx, guv, 1, len)
        idx += len
    end
    if !isnothing(controls.source.uq)
        gqv = vec(gq)
        len = length(gqv)
        copyto!(G, idx, gqv, 1, len)
        idx += len
    end
    if !isnothing(controls.massfrac.m)
        gmv = vec(gm)
        len = length(gmv)
        copyto!(G, idx, gmv, 1, len)
    end
    return nothing
end

"""
    prior_mass_fraction_cost(m, m₀, Qⁱₘ) -> Real

Calculate the quadratic cost associated with deviations of mass fractions `m` from a prior estimate `m₀`. This cost is weighted by the inverse covariance matrix `Qⁱₘ`, penalizing larger deviations from the prior. It is a core component of the objective function in an inverse problem.
"""
function prior_mass_fraction_cost(m::Vector,m₀::Vector,Qⁱₘ::Union{Diagonal, Symmetric})
    Δm = (m .- m₀)
    return (1 / length(m)) * Δm' * Qⁱₘ * Δm
end

function prior_mass_fraction_cost(m::NamedTuple,m₀::NamedTuple,Qⁱₘ::Union{Diagonal, Symmetric})
    return prior_mass_fraction_cost(vec(m), vec(m₀), Qⁱₘ)
end

"""
    gprior_mass_fraction_cost(m, m₀, Qₘ) -> Vector

Compute the gradient of the mass-fraction prior cost function. This gradient is essential for optimization algorithms that use first-derivative information to find the optimal state. It represents the sensitivity of the cost function to changes in the mass fractions.
"""
function gprior_mass_fraction_cost(m::Vector,m₀::Vector,Qₘ::Union{Diagonal, Symmetric})
    Δm = (m .- m₀)
    gmvec = (2 / length(m)) * Qₘ * Δm
    return gmvec
end

function gprior_mass_fraction_cost(m::NamedTuple,m₀::NamedTuple,Qₘ::Union{Diagonal, Symmetric})
    gmvec = gprior_mass_fraction_cost(vec(m), vec(m₀), Qₘ)
    return unvec(m₀, gmvec)
end
function gprior_mass_fraction_cost!(gm, m::NamedTuple,m₀::NamedTuple,Qₘ::Union{Diagonal, Symmetric})
    gmvec = gprior_mass_fraction_cost(vec(m), vec(m₀), Qₘ)
    adjustmassfraction!(gm, gmvec)
end

"""
    prior_boundary_cost(du, u₀, Qⁱᵤ) -> Real

Calculate the quadratic cost for perturbations `du` to the boundary conditions relative to a prior `u₀`. This term penalizes deviations from the prior boundary state, weighted by the inverse covariance `Qⁱᵤ`. It is used to regularize the inverse problem.
"""
function prior_boundary_cost(duvec::AbstractVector,
                             u₀::Union{BoundaryCondition{T}, NamedTuple},
                             Qⁱᵤ::Union{Diagonal, Symmetric}) where {T}
    return (1 / length(duvec)) * duvec' * Qⁱᵤ * duvec
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

Compute the gradient of the boundary condition prior cost. This provides the sensitivity of the cost function to changes in the boundary condition perturbations `du`. This gradient is a key input for gradient-based optimization of the model state.
"""
function gprior_boundary_cost(duvec::AbstractVector,
                              u₀::Union{BoundaryCondition{T}, NamedTuple},
                              Qⁱᵤ::Union{Diagonal, Symmetric}) where {T}
    return (2  / length(duvec)) * Qⁱᵤ * duvec
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

function gprior_boundary_cost!(gdu, du::NamedTuple, u₀::NamedTuple, Qᵤ::NamedTuple)
    for name in eachindex(du)
        gdu_i = gprior_boundary_cost(vec(du[name]), u₀[name], Qᵤ[name])
        adjustboundarycondition!(gdu[name], gdu_i; idx = 1, return_idx = false, r = 1.0)

    end

end

"""
    prior_source_cost(dq, q₀, Qⁱₛ) -> Real

Calculate the quadratic cost for perturbations `dq` to interior sources/sinks relative to a prior `q₀`. This term penalizes deviations from the prior source distribution, weighted by the inverse covariance `Qⁱₛ`. It helps regularize the solution of the inverse problem.
"""
function prior_source_cost(dqvec::Union{Vector, SubArray},
                           q₀::Union{Source{T}, NamedTuple},
                           Qⁱₛ::Union{Diagonal, Symmetric}) where {T}
    return (1 / length(dqvec)) * dqvec' * Qⁱₛ * dqvec
end

function prior_source_cost(
      dqvec::Vector,
      q₀::NamedTuple,
      Qₛ::NamedTuple{tracer_names,<:Tuple{Vararg{Union{Nothing, Symmetric,Diagonal}}}}
  ) where {tracer_names}
    J = 0
    lo = 1
    for name in tracer_names
        if !isnothing(q₀[name]) && !isnothing(get(Qₛ, name, nothing))
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

Compute the gradient of the interior source/sink prior cost. This provides the sensitivity of the cost function to changes in the source perturbations `dq`. This gradient is used by optimization algorithms to adjust the sources.
"""
function gprior_source_cost(dqvec::Union{Vector, SubArray},
                            q₀::Union{Source{T}, NamedTuple},
                            Qⁱₛ::Union{Diagonal, Symmetric}) where {T}
    gdqvec = (2 / length(dqvec))  * Qⁱₛ * dqvec
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

function gprior_source_cost!(gdq, dq::NamedTuple, q₀::NamedTuple, Qₛ::NamedTuple)
    for name in eachindex(dq)
        if !isnothing(Qₛ[name])
            gdq_i = gprior_source_cost(vec(dq[name]), q₀[name], Qₛ[name])
            adjustsource!(gdq[name], gdq_i)
        end
    end

end

"""
    model_observation_cost(n, Wⁱ) -> Real

Calculate the quadratic cost of the model-data misfit `n`, weighted by the inverse observation error covariance `Wⁱ`. This measures how poorly the model state reproduces the observations. Minimizing this term is the primary goal of the data assimilation.
"""
function model_observation_cost(n::Union{Vector, Field{T}},
                                Wⁱ::Union{Diagonal, Symmetric}) where {T}
    return (1 / length(vec(n))) * (n ⋅ (Wⁱ * n))
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

Compute the gradient of the model-observation misfit cost. This represents the sensitivity of the misfit cost to changes in the model state `n`. It is the initial step in the backward pass of the adjoint model.
"""
function gmodel_observation_cost(n::Union{Vector, Field{T}},
                                 Wⁱ::Union{Diagonal, Symmetric}) where {T}
    gnvec = (2 / length(vec(n))) * Wⁱ * n 
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
    model_data_misfit(c, c_obs, γ; ...) -> Field or Vector

Calculate the misfit (residual) between the modeled tracer field `c` and observations `c_obs`. This function handles both gridded and point observations, interpolating the model to observation locations if necessary. The result `n` is the difference `c - c_obs`.
"""
function model_data_misfit(c::Field{T}, c_obs::Field{T},
                                  γ::Grid; locs::Union{Nothing, Vector{G}}=nothing, wis = nothing) where {T, G}
    # n = (c - c_obs)
    n = c - c_obs
    return n
end

function model_data_misfit(c::Field{T}, c_obs::Vector,
                                  γ::Grid; locs::Union{Nothing, Vector{G}}=nothing, wis = nothing) where {T, G}
    # n = (c - c_obs)
    if !isnothing(wis)
        y = observe(c, wis, γ) 
    elseif !isnothing(locs)
        y = observe(c, locs, γ) 
    end
    
    n = y .- c_obs

    return n
end

function model_data_misfit(c::Field{T}, c_obs::Observations,
                                  γ::Grid; locs::Union{Nothing, Vector{G}}=nothing, wis = nothing) where {T, G}
    c_obs_val = c_obs.values
    n = model_data_misfit(c, c_obs_val, γ; locs=c_obs.locs) #need to specify wis here 
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
    gmodel_data_misfit(dn, c, c_obs, γ; ...) -> Field

Compute the adjoint of the model-data misfit calculation. This function propagates the sensitivities of the cost function with respect to the misfit `dn` back to sensitivities on the full tracer field `gc`. It is a crucial step in the adjoint model's backward pass.
"""
function gmodel_data_misfit(dn::Field{T}, c::Field{T}, c_obs::Field{T}, 
                                  γ::Grid; locs::Union{Nothing, NamedTuple}=nothing) where {T}
    gc = dn

    return gc

end

function gmodel_data_misfit(dn::Vector,
                                    c::Field{T}, c_obs::Vector, 
                                  γ::Grid; locs::Union{Nothing, NamedTuple}=nothing) where {T}
    gy = dn 
    gc = gobserve(gy,c,locs)  
    return gc

end

function gmodel_data_misfit(dn::Union{Vector, Field{T}},
                                    c::Field{T}, c_obs::Observations, 
                                  γ::Grid; locs::Union{Nothing, NamedTuple}=nothing) where {T}
    c_obs_val = c_obs.values
    if c_obs_val isa Field
        gc = dn
    elseif c_obs_val isa Vector
        gy = dn 
        gc = gobserve(gy,c,c_obs.locs; wis = c_obs.wis) 
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
    unconstrained_global_forward(controls, c_obs, γ; ...)

Execute the forward pass for the full cost function evaluation. This involves applying the current control parameters to solve for the steady-state tracer distribution and then calculating the total cost. It returns the total cost `J` and intermediate state variables needed for the corresponding backward pass.
"""
function unconstrained_global_forward(controls::Controls,
                                      c_obs, γ; locs = nothing)
    
    b = controls.boundary.b

    # Update all boundary fields (alyways independent)
    update_b!(controls.boundary)

    # Update all source fields (independent and coupled)
    update_q!(controls.source)
    q = controls.source.q

    A = watermassmatrix(controls.massfrac.m, γ, controls.massfrac.steps)
    # Alu = lu(A)
    prob = LinearProblem(A, Vector{Float64}(undef, size(A,1)))
    P = ilu(A; τ=0.01)
    cache = init(prob, KrylovJL_GMRES(); Pl=P)

    c = steadyinversion(cache, b, q, γ; c_obs=c_obs)

    # c = steadyinversion(Alu, b, q, γ)
    n = model_data_misfit(c, c_obs, γ; locs=locs)

    J_obs = model_observation_cost(n, c_obs) 
    J_source = prior_source_cost(controls.source.duq, controls.source.q₀, controls.source.Qₛ) 
    J_boundary = prior_boundary_cost(controls.boundary.dub, controls.boundary.u₀, controls.boundary.Qᵤ)
    J_massfrac = prior_mass_fraction_cost(controls.massfrac.m, controls.massfrac.m₀, controls.massfrac.Qₘ)
    # println("J_obs = ", J_obs)
    # println("J_source = ", J_source)
    # println("J_boundary = ", J_boundary)
    # println("J_massfrac = ", J_massfrac)
    J = J_obs + J_source + J_boundary + J_massfrac

    return (J=J, c=c, n=n, A=A)
end

"""
    unconstrained_global_backward!(state, controls, c_obs, γ; ...)

Execute the backward pass of the adjoint model to compute gradients of the cost function. This function takes the state from the forward pass and propagates sensitivities backward through the model equations. It calculates the gradients with respect to all control variables and stores them in the `controls` object.
"""
function unconstrained_global_backward!(state, controls::Controls, c_obs, γ; locs = nothing)
    c = state.c
    n = state.n
    A = state.A

    # Zero out all gradient buffers
    zero!(controls.boundary.gdub)
    zero!(controls.source.gduq)
    zero!(controls.source.gduq_cache)
    zero!(controls.massfrac.gm)

    # --- Prior Gradient Contributions ---
    gprior_boundary_cost!(controls.boundary.gdub, controls.boundary.dub, controls.boundary.u₀, controls.boundary.Qᵤ)
    gprior_source_cost!(controls.source.gduq, controls.source.duq, controls.source.q₀, controls.source.Qₛ)
    gprior_mass_fraction_cost!(controls.massfrac.gm, controls.massfrac.m, controls.massfrac.m₀, controls.massfrac.Qₘ)

    # --- Adjoint Model Pass ---
    # Start adjoint by differentiating the observation-cost wrt residuals.
    gn = gmodel_observation_cost(n, c_obs)
    gc = gmodel_data_misfit(gn, c, c_obs, γ; locs=locs)

    A_t = sparse(transpose(A))
    prob = LinearProblem(A_t, Vector{Float64}(undef, size(A_t,1)))
    P = ilu(A_t; τ=0.01)
    cache = init(prob, KrylovJL_GMRES(); Pl=P)

    # Calculate initial sensitivities wrt boundary and all sources (independent & dependent)
    gA_total = gsteadyinversion!(controls.boundary.gdub, controls.source.gduq_cache, 
                                 gc, c, A, cache, controls.boundary.b, controls.source.q, γ; c_obs=c_obs)

    # Apply chain rule for coupled sources
    update_gduq!(controls.source)

    # Propagate sensitivities back to mass fractions
    gwatermassmatrix!(controls.massfrac.gm, gA_total, controls.massfrac.m, γ, controls.massfrac.steps)
end



"""
    optim_fg_constrained_global_costfunction!(F, G, control_vector, ...)

This function provides the core interface to `Optim.jl` for cost function and gradient evaluation. It takes a flat `control_vector`, reconstructs the structured control parameters, runs the forward and backward passes, and returns the cost `F` and gradient `G`. It also handles transformations (e.g., softmax for mass fractions) required by the optimization.
"""
function optim_fg_constrained_global_costfunction!(F::Union{Nothing, Float64}, G::Union{Nothing, Vector{T}}, control_vector::Vector{T},
    controls::Controls, c_obs, γ; locs = nothing) where T

    unvec!(controls, control_vector) 
    #should make alpha an optional varibale, then check and apply transfomrationatins if needed
    α = 5. #an optional parameter to shrink gradients related to this transformation
    #controls.m is a real number vector 
    #transform m ∈ R to mass fraction that are non-negative and sum to 1
    softmax_massfractions!(controls.massfrac.m; α = α)  

    state = unconstrained_global_forward(controls, c_obs, γ; locs = locs)

    if G !== nothing
        unconstrained_global_backward!(state, controls, c_obs, γ; locs = locs)
        gx = gsoftmax_massfractions(controls.massfrac.gm, controls.massfrac.m; α = α)
        _write_gradient!(G, controls.boundary.gdub, controls.source.gduq, gx, controls)
        # GC.gc()
    end

    if F !== nothing
        return state.J
    end
end

function optim_fg_unconstrained_global_costfunction!(F::Union{Nothing, Float64}, G::Union{Nothing, Vector{T}}, control_vector::Vector{T},
    controls::Controls, c_obs, γ; locs = nothing, u₀_scale = nothing, q₀_scale = nothing ) where T

    unvec!(controls, control_vector) 

     # imposes a rescaling if desired by the user. in theory this should help improve the condition 
     # number of an approximated Hessian. However, modern solver (e.g., Ipopt) 
     # seem to deal with this quite well already.  
    if !isnothing(u₀_scale)
        rescale_parameter!(controls.boundary.ub, u₀_scale)
    end

    if !isnothing(q₀_scale)
        rescale_parameter!(controls.source.uq, q₀_scale)
    end
    
    state = unconstrained_global_forward(controls, c_obs, γ; locs = locs)

    if G !== nothing
        unconstrained_global_backward!(state, controls, c_obs, γ; locs = locs)

        if !isnothing(u₀_scale)
            rescale_parameter!(controls.boundary.gdub, u₀_scale)
        end

        if !isnothing(q₀_scale)
            rescale_parameter!(controls.source.gduq, q₀_scale)
        end       

        _write_gradient!(G, controls.boundary.gdub, controls.source.gduq, controls.massfrac.gm, controls)
    end

    if F !== nothing
        return state.J
    end
end
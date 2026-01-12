
"""
    joint_global_cost!(F, G, control_vector, ...)

    Interface for `Optim.jl` for a 'joint' inversion on a 'global' grid.

    A 'joint' inversion simultaneously optimizes all control parameters (boundary 
    conditions, sources, and mass fractions). A 'global' grid refers to the 
    entire world ocean. This function assumes an unconstrained optimization problem.

# Arguments
- `F::Union{Nothing,Float64}`: cost accumulator; when `nothing`, only the gradient is computed.
- `G::Union{Nothing,Vector}`: gradient buffer populated in-place when provided.
- `control_vector::Vector`: flat control vector passed to the optimizer.
- `controls::Controls`: structured control container that mirrors `control_vector`.
- `c_obs`: observations for the misfit term.
- `γ::Grid`: model grid.
- `locs`: optional observation locations for point data.
- `u₀_scale`: optional scaling for boundary-condition controls and gradients.
- `q₀_scale`: optional scaling for source controls and gradients.
# Output
- `Float64` or `nothing`: returns the cost when `F` is not `nothing`.
"""
function joint_global_cost!(F::Union{Nothing, Float64}, 
    G::Union{Nothing, Vector{T}}, control_vector::Vector{T},
    controls::Controls, c_obs, γ; locs = nothing, 
    u₀_scale = nothing, q₀_scale = nothing, 
    debug = false) where T

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
    
    state = forward_joint_global(controls, c_obs, γ; locs = locs, debug = debug)

    if G !== nothing
        backward_joint_global!(state, controls, c_obs, γ; locs = locs)

        if !isnothing(u₀_scale)
            rescale_parameter!(controls.boundary.gdub, u₀_scale)
        end

        if !isnothing(q₀_scale)
            rescale_parameter!(controls.source.gduq, q₀_scale)
        end       

        write_gradient!(G, controls.boundary.gdub, controls.source.gduq, controls.massfrac.gm, controls)
    end

    if F !== nothing
        return state.J
    end
end

"""
    forward_joint_global(controls, c_obs, γ; ...)

    Forward pass for the 'joint global' problem.

    Applies all controls (boundary conditions, sources, mass fractions), 
    solves for the tracer field, and accumulates all cost components.
# Arguments
- `controls::Controls`: current control parameters for boundary, source, and mass fractions.
- `c_obs`: observations (gridded or point) used in the misfit term.
- `γ::Grid`: model grid.
- `locs`: optional observation locations when working with point data.
# Output
- `NamedTuple`: `(J, c, n, A)` containing total cost, modeled tracer, residuals, and system matrix.
"""
function forward_joint_global(controls::Controls,
                                      c_obs, γ; locs = nothing, debug = false)
    
    b = controls.boundary.b

    # Update all boundary fields (alyways independent)
    update_b!(controls.boundary)

    # Update all source fields (independent and coupled)
    update_q!(controls.source)
    q = controls.source.q

    A = watermassmatrix(controls.massfrac.m, γ, controls.massfrac.steps)
    Alu = lu(A)
    # prob = LinearProblem(A, Vector{Float64}(undef, size(A,1)))
    # P = ilu(A; τ=0.01)
    # cache = init(prob, KrylovJL_GMRES(); Pl=P)
    # c = steadyinversion(cache, b, q, γ; c_obs=c_obs)

    c = steadyinversion(Alu, b, q, γ)
    n = model_data_misfit(c, c_obs, γ; locs=locs)

    J_obs = model_observation_cost(n, c_obs) 
    J_source = prior_source_cost(controls.source.duq, controls.source.q₀, controls.source.Qₛ) 
    J_boundary = prior_boundary_cost(controls.boundary.dub, controls.boundary.u₀, controls.boundary.Qᵤ)
    J_massfrac = prior_mass_fraction_cost(controls.massfrac.m, controls.massfrac.m₀, controls.massfrac.Qₘ)
    J = J_obs + J_source + J_boundary + J_massfrac

    if debug 
        println("J_obs = ", J_obs)
        println("J_source = ", J_source)
        println("J_boundary = ", J_boundary)
        println("J_massfrac = ", J_massfrac)
    end
    
    return (J=J, c=c, n=n, A=A, Alu = Alu)
end

"""
    backward_joint_global!(state, controls, c_obs, γ; ...)

    Backward pass for the 'joint global' problem.

    Computes gradients for all control blocks (boundary conditions, sources,
    and mass fractions) via the adjoint model.
# Arguments
- `state`: output of `forward_joint_global` containing `c`, `n`, and `A`.
- `controls::Controls`: control container holding gradient buffers to be filled.
- `c_obs`: observations used to weight the misfit.
- `γ::Grid`: model grid.
- `locs`: optional observation locations for point data.
# Output
- `nothing`: populates gradient fields stored inside `controls`.
"""
function backward_joint_global!(state, controls::Controls, c_obs, γ; locs = nothing)
    c = state.c
    n = state.n
    A = state.A
    Alu = state.Alu

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

    # A_t = sparse(transpose(A))
    # prob = LinearProblem(A_t, Vector{Float64}(undef, size(A_t,1)))
    # P = ilu(A_t; τ=0.01)
    # cache = init(prob, KrylovJL_GMRES(); Pl=P)

    # # Calculate initial sensitivities wrt boundary and all sources (independent & dependent)
    # gA_total = gsteadyinversion!(controls.boundary.gdub, controls.source.gduq_cache, 
    #                              gc, c, A, cache, controls.boundary.b, controls.source.q, γ; c_obs=c_obs)
    gA_total = gsteadyinversion!(controls.boundary.gdub, controls.source.gduq_cache, 
                                 gc, c, A, Alu, controls.boundary.b, controls.source.q, γ; c_obs=c_obs)
    # Apply chain rule for coupled sources
    update_gduq!(controls.source)

    # Propagate sensitivities back to mass fractions
    gwatermassmatrix!(controls.massfrac.gm, gA_total, controls.massfrac.m, γ, controls.massfrac.steps)
end




"""
    prior_mass_fraction_cost(m, m₀, Qⁱₘ) -> Real
    Quadratic prior penalty for the mass fractions.
# Arguments
- `m::Union{Vector,NamedTuple}`: current mass fractions (flattened or structured).
- `m₀::Union{Vector,NamedTuple}`: prior mass fractions in the same layout as `m`.
- `Qⁱₘ::Union{Diagonal,Symmetric}`: inverse covariance weighting for the prior.
# Output
- `Real`: mean-square cost for deviations from the prior mass fractions.
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
    Gradient of the mass-fraction prior penalty.
# Arguments
- `m::Union{Vector,NamedTuple}`: current mass fractions.
- `m₀::Union{Vector,NamedTuple}`: prior mass fractions.
- `Qₘ::Union{Diagonal,Symmetric}`: covariance (not inverse) used to scale the gradient.
# Output
- `Vector` or `NamedTuple`: gradient with respect to `m`, matching the layout of the input.
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
    Quadratic prior penalty for boundary-condition perturbations.
# Arguments
- `du::Union{AbstractVector,NamedTuple}`: perturbations to boundary conditions, flattened or structured.
- `u₀::Union{BoundaryCondition,NamedTuple}`: prior boundary state.
- `Qⁱᵤ::Union{Diagonal,Symmetric,NamedTuple}`: inverse covariance weighting for each boundary control.
# Output
- `Real`: mean-square cost for deviating boundary controls from their priors.
"""
function prior_boundary_cost(duvec::AbstractVector,
                             u₀::Union{BoundaryCondition{T}, NamedTuple},
                             Qⁱᵤ::Union{Diagonal, Symmetric}) where {T}
    return (1 / length(duvec)) * (duvec' * Qⁱᵤ * duvec)
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
    println("\nStarting Boundary Cost")
    for (i, name) in enumerate(tracer_names)
        J1 = prior_boundary_cost(vec(du[name]), u₀[name], Qᵤ[name])
        println("$name : $J1")
        J += J1
    end
    println("Ending Boundary Cost")

    return J
end

"""
    gprior_boundary_cost(du, u₀, Qⁱᵤ) -> Vector
    Gradient of the boundary prior penalty.
# Arguments
- `du::Union{AbstractVector,NamedTuple}`: boundary-condition perturbations.
- `u₀::Union{BoundaryCondition,NamedTuple}`: prior boundary state (used for sizing).
- `Qⁱᵤ::Union{Diagonal,Symmetric,NamedTuple}`: inverse covariance weighting.
# Output
- `Vector` or `NamedTuple`: gradient with respect to `du`, matching the layout of the input.
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
    Quadratic prior penalty for interior sources and sinks.
# Arguments
- `dq::Union{Vector,SubArray,NamedTuple}`: perturbations to interior sources.
- `q₀::Union{Source,NamedTuple}`: prior source field (may contain `nothing` for inactive tracers).
- `Qⁱₛ::Union{Diagonal,Symmetric,NamedTuple}`: inverse covariance weighting for source controls.
# Output
- `Real`: mean-square cost for deviating interior sources from their priors.
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
    Gradient of the interior source prior penalty.
# Arguments
- `dq::Union{Vector,SubArray,NamedTuple}`: interior-source perturbations.
- `q₀::Union{Source,NamedTuple}`: prior source field.
- `Qⁱₛ::Union{Diagonal,Symmetric,NamedTuple}`: inverse covariance weighting.
# Output
- `Vector` or `NamedTuple`: gradient with respect to `dq`, matching the input layout.
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
    Mean-square misfit between model results and observations.
# Arguments
- `n::Union{Vector,Field,NamedTuple}`: model-data residuals.
- `Wⁱ::Union{Diagonal,Symmetric,NamedTuple}`: inverse observation-error covariance.
# Output
- `Real`: cost term measuring the observation misfit.
"""
function model_observation_cost(n::Union{Vector, Field{T}},
                                Wⁱ::Union{Diagonal, Symmetric}) where {T}
    if n isa Field
        # Create a copy and zero out non-interior points
        n_mod = deepcopy(n)
        boundary_mask = n_mod.γ.wet .& .!n_mod.γ.interior
        n_mod.tracer[boundary_mask] .*= 0.0
        norm_factor = sum(n_mod.γ.interior)
        
        return (1 / norm_factor) * n_mod ⋅ (Wⁱ * n_mod)
    else
        # Keep original behavior for vector observations
        return (1 / length(vec(n))) * (n ⋅ (Wⁱ * n))
    end
end

function model_observation_cost(n::NamedTuple,c_obs::NamedTuple)
    tracer_names = keys(n)
    J = 0
    println("Starting ModelObs Cost")
    for (i, name) in enumerate(tracer_names)
        J1 = model_observation_cost(n[name], c_obs[name].W)
        println("$name : $J1")
        J += J1
    end
    println("Ending ModelObs Cost")
    return J
end

"""
    gmodel_observation_cost(n, Wⁱ) -> Vector
    Gradient of the observation misfit cost with respect to residuals.
# Arguments
- `n::Union{Vector,Field,NamedTuple}`: model-data residuals.
- `Wⁱ::Union{Diagonal,Symmetric,NamedTuple}`: inverse observation-error covariance.
# Output
- `Vector` or `NamedTuple`: gradient with respect to `n`, matching the input layout.
"""
function gmodel_observation_cost(n::Union{Vector, Field{T}},
                                 Wⁱ::Union{Diagonal, Symmetric}) where {T}
    if n isa Field
        # Create a copy and zero out non-interior points
        n_mod = deepcopy(n)
        boundary_mask = n_mod.γ.wet .& .!n_mod.γ.interior
        n_mod.tracer[boundary_mask] .*= 0.0
        norm_factor = sum(n_mod.γ.interior)
        return (2 / norm_factor) * (Wⁱ * n_mod)

    else
        return (2 / length(vec(n))) * (Wⁱ * n)
    end
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
    Compute residuals between modeled and observed tracer fields.
# Arguments
- `c::Union{Field,NamedTuple}`: modeled tracer field.
- `c_obs::Union{Field,Vector,Observations,NamedTuple}`: observations, potentially with locations and weights.
- `γ::Grid`: model grid used to evaluate or interpolate the field.
- `locs`: optional observation locations; when omitted, uses locations stored in `c_obs`.
- `wis`: optional interpolation weights for mapping from the grid to observations.
# Output
- `Field`, `Vector`, or `NamedTuple`: residuals `c - c_obs` in the same layout as `c`.
"""
function model_data_misfit(c::Field{T}, c_obs::Field{T},
                                  γ::Grid; locs::Union{Nothing, Vector{G}}=nothing, wis = nothing) where {T, G}
    # n = (c - c_obs)
    n = c - c_obs; 

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
    Adjoint of the model-data misfit, mapping residual gradients to tracer gradients.
# Arguments
- `dn::Union{Field,Vector,NamedTuple}`: sensitivities with respect to the residuals.
- `c::Union{Field,NamedTuple}`: modeled tracer field used for sizing/interpolation.
- `c_obs::Union{Field,Vector,Observations,NamedTuple}`: observations or their metadata.
- `γ::Grid`: model grid.
- `locs`: optional observation locations; when omitted, uses those stored in `c_obs`.
# Output
- `Field` or `NamedTuple`: gradient with respect to the modeled tracer field `c`.
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


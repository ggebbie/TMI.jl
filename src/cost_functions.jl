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
function gprior_mass_fraction_cost!(gm, m::NamedTuple,m₀::NamedTuple,Qₘ::Union{Diagonal, Symmetric})
    gmvec = gprior_mass_fraction_cost(vec(m), vec(m₀), Qₘ)
    adjustmassfraction!(gm, gmvec)
end


"""
    prior_boundary_cost(du, u₀, Qⁱᵤ) -> Real

Quadratic penalty on boundary-condition perturbations: `du' Qⁱᵤ du`.
"""
function prior_boundary_cost(duvec::AbstractVector,
                             u₀::Union{BoundaryCondition{T}, NamedTuple},
                             Qⁱᵤ::Union{Diagonal, Symmetric}) where {T}
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
                              u₀::Union{BoundaryCondition{T}, NamedTuple},
                              Qⁱᵤ::Union{Diagonal, Symmetric}) where {T}
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

function gprior_boundary_cost!(gdu, du::NamedTuple, u₀::NamedTuple, Qᵤ::NamedTuple)
    for name in eachindex(du)
        gdu_i = gprior_boundary_cost(vec(du[name]), u₀[name], Qᵤ[name])
        adjustboundarycondition!(gdu[name], gdu_i; idx = 1, return_idx = false, r = 1.0)

    end

end


"""
    prior_source_cost(dq, q₀, Qⁱₛ) -> Real

Quadratic penalty on source perturbations: `dq' Qⁱₛ dq`.
"""
function prior_source_cost(dqvec::Union{Vector, SubArray},
                           q₀::Union{Source{T}, NamedTuple},
                           Qⁱₛ::Union{Diagonal, Symmetric}) where {T}
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
function gprior_source_cost(dqvec::Union{Vector, SubArray},
                            q₀::Union{Source{T}, NamedTuple},
                            Qⁱₛ::Union{Diagonal, Symmetric}) where {T}
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

Weighted observation misfit: `n' Wⁱ n`.
"""
function model_observation_cost(n::Union{Vector, Field{T}},
                                Wⁱ::Union{Diagonal, Symmetric}) where {T}
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
function gmodel_observation_cost(n::Union{Vector, Field{T}},
                                 Wⁱ::Union{Diagonal, Symmetric}) where {T}
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
    gmodel_data_misfit(dn, c, c_obs, γ; locs=nothing) -> gradient

Adjoint of `model_data_misfit`: propagate residual sensitivities `dn` back to the state.
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
    unconstrained_global_forward(ub, uq, m, controls, c_obs, γ; locs=nothing)

Forward pass for the unconstrained cost: applies control updates, assembles matrices,
solves the steady state, and returns intermediates needed for the backward pass.
"""
function unconstrained_global_forward(controls::ControlParameters,
                                      c_obs, γ; locs = nothing)
    b = controls.boundary.b # working boundary fields
    q = controls.source.q

    du_names = keys(controls.boundary.ub)
    uq_base_keys = _uq_base_keys(controls.source.uq)

    @inbounds for key in keys(b)
        setboundarycondition!(controls.boundary.b[key], controls.boundary.u₀[key]) #du = ub - u₀
        if key ∈ du_names
            setboundarycondition!(controls.boundary.dub[key], controls.boundary.ub[key]) #du = ub - u₀
            adjustboundarycondition!(controls.boundary.dub[key], controls.boundary.u₀[key]; r = -1.0) #du = ub - u₀
            adjustboundarycondition!(b[key], controls.boundary.dub[key]; r = 1.0) #b = b + du
        end
    end

    @inbounds for key in keys(q)
        if isnothing(q[key])
            continue
        end
        q0_key = isnothing(controls.source.q₀) ? nothing : controls.source.q₀[key]
        if isnothing(q0_key)
            zero!(controls.source.q[key])                 # child with no prior: start from zero
        else
            replacesource!(controls.source.q[key], q0_key) # dq = uq - q₀
        end
        if key ∈ uq_base_keys #mirrors existing behavior; keys(du) drives source updates
            replacesource!(controls.source.duq[key], controls.source.uq[key]) #du = ub - u₀
            if !isnothing(q0_key)
                adjustsource!(controls.source.duq[key], q0_key; r = -1.0) #dq = uq - q₀
            end
            adjustsource!(q[key], controls.source.duq[key]; r = 1.0)  #q = q + dq
        end
    end

    # Apply dependent source links: dq_child = scale * dq_parent
    for (child, link) in pairs(_uq_links(controls.source.uq))
        parent = link.dependson
        scale = link.scale
        zero!(controls.source.duq[child])
        adjustsource!(controls.source.duq[child], controls.source.duq[parent]; r = scale)
        adjustsource!(q[child], controls.source.duq[child]; r = 1.0)
    end

    A = watermassmatrix(controls.massfrac.m, γ, controls.massfrac.mass_fraction_steps)
    # Alu = lu(A) #avoiding the LU decomposition for now. 

    prob = LinearProblem(A,  Vector{Float64}(undef, size(A,1)))
    P = ilu(A; τ=0.01)
    cache = init(prob, KrylovJL_GMRES(); Pl=P)

    # c = steadyinversion(Alu, b, q, γ)
    c =  steadyinversion(cache, b, q, γ)
    n = model_data_misfit(c, c_obs, γ; locs=locs)

    J = model_observation_cost(n,c_obs) 
    J += prior_source_cost(controls.source.duq, controls.source.q₀, controls.source.Qₛ) 
    J += prior_boundary_cost(controls.boundary.dub, controls.boundary.u₀, controls.boundary.Qᵤ)
    J += prior_mass_fraction_cost(controls.massfrac.m,controls.massfrac.m₀,controls.massfrac.Qₘ)

    return (
        J = J,
        c = c,
        n = n,
        A = A,
    )
end

"""
    unconstrained_global_backward(state, controls, c_obs, γ; locs=nothing)

Backward pass given intermediates from `unconstrained_global_forward`, returning gradients
with respect to (ub, uq, m) in the same structured containers as their inputs.
"""
function unconstrained_global_backward!(state, controls::ControlParameters, c_obs, γ; locs = nothing)
    c = state.c
    n = state.n
    A = state.A

    zero!(controls.boundary.gdub)
    zero!(controls.source.gduq)
    zero!(controls.massfrac.gm)

    gprior_boundary_cost!(controls.boundary.gdub, controls.boundary.dub,controls.boundary.u₀, controls.boundary.Qᵤ)
    gprior_source_cost!(controls.source.gduq, controls.source.duq,controls.source.q₀, controls.source.Qₛ)
    gprior_mass_fraction_cost!(controls.massfrac.gm, controls.massfrac.m, controls.massfrac.m₀, controls.massfrac.Qₘ)

    # Start adjoint by differentiating the observation-cost wrt residuals.
    gn = gmodel_observation_cost(n, c_obs)
    gc = gmodel_data_misfit(gn, c, c_obs, γ; locs=locs)

    A_t = sparse(transpose(A))
    prob = LinearProblem(A_t,Vector{Float64}(undef, size(A_t,1)))
    P = ilu(A_t; τ=0.01)
    cache = init(prob, KrylovJL_GMRES(); Pl=P)

    gA_total = gsteadyinversion!(controls.boundary.gdub, controls.source.gduq, 
                                 gc, c, A, cache, controls.boundary.b, controls.source.q, γ)

    gwatermassmatrix!(controls.massfrac.gm, gA_total, controls.massfrac.m, γ, controls.massfrac.mass_fraction_steps)

    # Accumulate dependent-source gradients back to their bases.
    for (child, link) in pairs(_uq_links(controls.source.uq))
        parent = link.dependson
        scale = link.scale
        adjustsource!(controls.source.gduq[parent], controls.source.gduq[child]; r = scale)
    end

end

@inline function _write_gradient!(G, gu, gq, gm, controls::ControlParameters)
    idx = 1
    if !isnothing(controls.boundary.ub)
        guv = vec(gu)
        len = length(guv)
        copyto!(G, idx, guv, 1, len)
        idx += len
    end
    if !isnothing(controls.source.uq)
        base_keys = _uq_base_keys(controls.source.uq)
        gqv = _vec_base_sources(gq, base_keys)
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

function optim_fg_constrained_global_costfunction!(F::Union{Nothing, Float64}, G::Union{Nothing, Vector{T}}, control_vector::Vector{T},
    controls::ControlParameters, c_obs, γ; locs = nothing) where T

    unvec!(controls, control_vector) 
    #should make alpha an optional varibale, then check and apply transfomrationatins if needed
    α = 5. #an optional parameter to shrink gradients related to this transformation
    # controls.massfrac.m is a real-number vector; transform to non-negative fractions summing to 1
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

###### OLD REFERENCE. More anaytical, but really slow #####
# """
#     unconstrained_global_costfunction(ub, uq, m, controls, c_obs, γ; locs=nothing, return_gradients=false)

# Objective for the unconstrained inversion: update `b`/`q` with controls, assemble `A`,
# solve the steady state, add observation misfit and priors, and optionally return adjoint
# gradients w.r.t. `du`, `dq`, and `m` when `return_gradients=true`.
# """
# function unconstrained_global_costfunction(ub, uq, m, controls::ControlParameters, 
#                                            c_obs, γ; locs = nothing, return_gradients = false)
#     # Forward pass: compute state at current control vector
#     b = deepcopy(controls.u₀) #deep-copying to 
#     q = deepcopy(controls.q₀)
#     du = deepcopy(ub)
#     dq = deepcopy(uq)

#     # Apply control updates: du,dq store perturbations relative to priors.
#     du_names = keys(ub)
#     dq_names = keys(uq)

#     for key in keys(b)
#         if key ∈ du_names
#             adjustboundarycondition!(du[key], controls.u₀[key]; r = -1.0) #du = ub - u₀
#             adjustboundarycondition!(b[key], du[key]; r = 1.0) #du = ub - u₀
#         end
#     end

#     for key in keys(q)
#         if !isnothing(q[key])
#             if key ∈ du_names
#                 adjustsource!(dq[key], controls.q₀[key]; r = -1.0) #dq = uq - q₀
#                 adjustsource!(q[key], dq[key]; r = 1.0)  #uq = uq + q₀
#             end
#         end
#     end

#     A = watermassmatrix(m, γ)
#     Alu = lu(A)

#     # Solve for each tracer
#     c = steadyinversion(Alu,b, q, γ)
#     n = model_data_misfit(c, c_obs, γ; locs=locs)
#     J = model_observation_cost(n,c_obs) + 
#         prior_source_cost(dq, controls.q₀, controls.Qₛ) + 
#         prior_boundary_cost(du, controls.u₀, controls.Qᵤ) + 
#         prior_mass_fraction_cost(m,controls.m₀,controls.Qₘ)
        
#     if !return_gradients 
#         return J, nothing, nothing, nothing #just return cost
#     else
#         # Propagate gradients through adjoint operations
#         gdu = zero(du)
#         gdq = zero(dq)
#         gm = similar(controls.m₀)

#         # Prior contributions (no model dynamics).
#         gdu_1 = gprior_boundary_cost(du,controls.u₀, controls.Qᵤ)
#         gdq_1 = gprior_source_cost(dq,controls.q₀, controls.Qₛ)
#         gm_1 = gprior_mass_fraction_cost(m, controls.m₀, controls.Qₘ)

#         # Start adjoint by differentiating the observation-cost wrt residuals.
#         gn = gmodel_observation_cost(n, c_obs)
#         gc = gmodel_data_misfit(gn, c, c_obs, γ; locs=locs)
#         gdu_2, gdq_2, gA_total = gsteadyinversion(gc, c, A, Alu, b, q)
#         gm_2 = gwatermassmatrix(gA_total, m, γ)

#         # Accumulate adjoint contributions into control gradients.
#         for key in du_names
#             gadjustboundarycondition!(gdu[key], gdu_1[key])
#             gadjustboundarycondition!(gdu[key], gdu_2[key])
#         end

#         for key in dq_names
#             if !isnothing(q[key])
#                 gadjustsource!(gdq[key], gdq_1[key], controls.q₀[key])
#                 gadjustsource!(gdq[key], gdq_2[key], controls.q₀[key])
#             end
#         end

#         for key in keys(gm)
#             gm[key].fraction .= gm_1[key].fraction .+ gm_2[key].fraction
#         end

#         #since du = u - u₀, ∂J/∂du = ∂J/∂u. Same for dq if q is not on logscale 
#         return J, gdu, gdq, gm
#     end
# end

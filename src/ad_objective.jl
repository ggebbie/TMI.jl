struct ADTemplates{B,S,M,D,G}
    boundary_prior::B
    source_prior::S
    massfrac_steps::M
    dependencies::D
    γ::G
end

ADTemplates(controls::Controls) = ADTemplates(
    controls.boundary.u₀,
    controls.source.q₀,
    controls.massfrac.steps,
    controls.source.dependencies,
    controls.γ,
)

_subset_namedtuple(nt::NamedTuple, names) = NamedTuple{names}(map(name -> getproperty(nt, name), names))

function _boundary_from_controls(ub::NamedTuple, u₀::NamedTuple)
    tracer_names = keys(u₀)
    values = map(tracer_names) do name
        haskey(ub, name) ? ub[name] : u₀[name]
    end
    return NamedTuple{tracer_names}(values)
end

function _source_from_controls(uq::NamedTuple, q₀::NamedTuple, dependencies::NamedTuple)
    independent = keys(uq)
    all_names = keys(q₀)
    q_values = map(all_names) do name
        if name in independent
            uq[name]
        elseif haskey(dependencies, name)
            indep_name, alpha = dependencies[name]
            alpha * uq[indep_name]
        else
            q₀[name]
        end
    end
    return NamedTuple{all_names}(q_values)
end

function _source_control_gradient(gq_all::NamedTuple, uq::NamedTuple, dependencies::NamedTuple)
    independent = keys(uq)
    grad_values = map(independent) do name
        grad = deepcopy(gq_all[name])
        for dep_name in keys(dependencies)
            indep_name, alpha = dependencies[dep_name]
            if indep_name == name
                adjustsource!(grad, gq_all[dep_name]; r = alpha)
            end
        end
        grad
    end
    return NamedTuple{independent}(grad_values)
end

function _scaled_control_vector(
    x::AbstractVector,
    controls::Controls;
    u₀_scale = nothing,
    q₀_scale = nothing,
)
    idx = 1
    parts = Vector{eltype(x)}[]

    if !isnothing(controls.boundary.ub)
        len = controls.ub_length
        x_part = x[idx:idx+len-1]
        if isnothing(u₀_scale)
            push!(parts, x_part)
        else
            scale_vec = vec(_subset_namedtuple(u₀_scale, keys(controls.boundary.ub)))
            push!(parts, x_part .* scale_vec)
        end
        idx += len
    end

    if !isnothing(controls.source.uq)
        len = controls.uq_length
        x_part = x[idx:idx+len-1]
        if isnothing(q₀_scale)
            push!(parts, x_part)
        else
            # Match current behavior in `joint_global_cost!`, where source scaling is
            # routed through `rescale_parameter!` and currently has no effect on `Source`.
            push!(parts, x_part)
        end
        idx += len
    end

    if !isnothing(controls.massfrac.m)
        push!(parts, x[idx:end])
    end

    return isempty(parts) ? similar(x, 0) : vcat(parts...)
end

function _boundary_prior_cost_ad(ub::NamedTuple, controls::Controls)
    u₀_active = _subset_namedtuple(controls.boundary.u₀, keys(ub))
    du = vec(ub) .- vec(u₀_active)
    return prior_boundary_cost(du, controls.boundary.u₀, controls.boundary.Qᵤ)
end

function _source_prior_cost_ad(uq::NamedTuple, controls::Controls)
    q₀_active = _subset_namedtuple(controls.source.q₀, keys(uq))
    dq = vec(uq) .- vec(q₀_active)
    return prior_source_cost(dq, controls.source.q₀, controls.source.Qₛ)
end

function _massfrac_prior_cost_ad(m::NamedTuple, controls::Controls)
    return prior_mass_fraction_cost(vec(m), vec(controls.massfrac.m₀), controls.massfrac.Qₘ)
end

function _model_observation_cost_ad(n::Field{T}, Wⁱ::Union{Diagonal, Symmetric}) where {T}
    nvec = vec(n)
    mask = n.γ.interior[n.γ.wet]
    masked = nvec .* mask
    norm_factor = sum(n.γ.interior)
    return (1 / norm_factor) * dot(masked, Wⁱ * masked)
end

function _model_observation_cost_ad(n::AbstractVector, Wⁱ::Union{Diagonal, Symmetric})
    return (1 / length(n)) * dot(n, Wⁱ * n)
end

function _model_observation_cost_ad(n::NamedTuple, c_obs::NamedTuple)
    tracer_names = keys(n)
    total = zero(Float64)
    for name in tracer_names
        total += _model_observation_cost_ad(n[name], c_obs[name].W)
    end
    return total
end

function _field_like(template::Field, tracer)
    return Field(tracer, template.γ, template.name, template.longname, template.units)
end

function _solve_from_controls(
    ub::NamedTuple,
    uq::NamedTuple,
    m::NamedTuple,
    templates::ADTemplates,
)
    b = _boundary_from_controls(ub, templates.boundary_prior)
    q = _source_from_controls(uq, templates.source_prior, templates.dependencies)
    A = watermassmatrix(m, templates.γ, templates.massfrac_steps)
    Alu = lu(A)
    c = steadyinversion(Alu, b, q, templates.γ)
    return c, (b = b, q = q, A = A, Alu = Alu, c = c)
end

function ChainRulesCore.rrule(
    ::typeof(_solve_from_controls),
    ub::NamedTuple,
    uq::NamedTuple,
    m::NamedTuple,
    templates::ADTemplates,
)
    c, cache = _solve_from_controls(ub, uq, m, templates)

    function solve_pullback(Δc)
        if Δc isa AbstractZero
            return NoTangent(), zero(ub), zero(uq), zero(m), NoTangent()
        end

        gc = if Δc isa NamedTuple
            Δc
        else
            ProjectTo(c)(Δc)
        end

        gdub = zero(ub)
        gq_all = zero(cache.q)
        gA_total = gsteadyinversion!(
            gdub,
            gq_all,
            gc,
            cache.c,
            cache.A,
            cache.Alu,
            cache.b,
            cache.q,
            templates.γ;
            c_obs = NamedTuple(),
        )

        gduq = _source_control_gradient(gq_all, uq, templates.dependencies)
        gm = gwatermassmatrix(gA_total, m, templates.γ, templates.massfrac_steps)

        return NoTangent(), gdub, gduq, gm, NoTangent()
    end

    return c, solve_pullback
end

function joint_global_objective_ad(
    control_vector::AbstractVector,
    controls::Controls,
    c_obs,
    γ;
    locs = nothing,
    u₀_scale = nothing,
    q₀_scale = nothing,
    debug = false,
)
    x_scaled = _scaled_control_vector(
        control_vector,
        controls;
        u₀_scale = u₀_scale,
        q₀_scale = q₀_scale,
    )
    ub, uq, m = unvec(controls, collect(x_scaled))
    templates = ADTemplates(controls)
    c = _solve_from_controls(ub, uq, m, templates)[1]
    n = model_data_misfit(c, c_obs, γ; locs = locs)

    J_obs = _model_observation_cost_ad(n, c_obs)
    J_source = _source_prior_cost_ad(uq, controls)
    J_boundary = _boundary_prior_cost_ad(ub, controls)
    J_massfrac = _massfrac_prior_cost_ad(m, controls)

    if debug
        println("J_obs = ", J_obs)
        println("J_source = ", J_source)
        println("J_boundary = ", J_boundary)
        println("J_massfrac = ", J_massfrac)
    end

    return J_obs + J_source + J_boundary + J_massfrac
end

function joint_global_cost_ad!(
    F::Union{Nothing, Float64},
    G::Union{Nothing, Vector{T}},
    control_vector::Vector{T},
    controls::Controls,
    c_obs,
    γ;
    locs = nothing,
    u₀_scale = nothing,
    q₀_scale = nothing,
    debug = false,
) where {T}
    f = x -> joint_global_objective_ad(
        x,
        controls,
        c_obs,
        γ;
        locs = locs,
        u₀_scale = u₀_scale,
        q₀_scale = q₀_scale,
        debug = debug,
    )

    if G !== nothing
        result = Zygote.withgradient(f, control_vector)
        G .= result.grad[1]
        return F === nothing ? nothing : result.val
    end

    return F === nothing ? nothing : f(control_vector)
end

function check_shared_references(nt::NamedTuple, name::String)
    vals = filter(!isnothing, collect(values(nt)))
    n = length(vals)
    for i in 1:n, j in (i+1):n
        if vals[i] === vals[j]
            key_names = collect(keys(nt))
            shared_keys = String[]
            for k in 1:length(nt)
                if !isnothing(nt[k]) && (nt[k] === vals[i])
                    push!(shared_keys, string(key_names[k]))
                end
            end
            msg = """Shared reference detected in $name: fields $(join(shared_keys, ", ")) point to the same object in memory. Use deepcopy to create independent copies."""
            error(msg)
        end
    end
    return nothing
end

check_shared_references(x, ::String) = nothing # Non-NamedTuple fallback

check_bounds_keys(::Nothing, ::Any, ::String) = nothing
function check_bounds_keys(bounds::NamedTuple, controls::NamedTuple, name::String)
    keys(bounds) == keys(controls) || error("$name keys must match control keys")
    return nothing
end

struct ControlParameters{BC<:BoundaryControls,SC<:SourceControls,MC<:MassFractionControls,V,UBL,QL,ML,G<:Grid}
    boundary::BC
    source::SC
    massfrac::MC
    vector::V
    ub_length::UBL
    uq_length::QL
    m_length::ML
    γ::G
end

function ControlParameters(; γ::Grid,
                            ub = nothing, uq = nothing, m = nothing,
                            u₀ = nothing, q₀ = nothing, m₀ = nothing,
                            Qᵤ = nothing, Qₛ = nothing, Qₘ = nothing,
                            ub_bounds = nothing, uq_bounds = nothing, m_bounds = nothing,
                            cache_vector = false,
                            cache_mass_fraction_steps = true,
                            cache_A = true)
    # Ensure control entries and priors are independent objects.
    for (val, name) in ((ub, "ub"), (uq, "uq"), (m, "m"),
                        (u₀, "u₀"), (q₀, "q₀"), (m₀, "m₀"),
                        (Qᵤ, "Qᵤ"), (Qₛ, "Qₛ"))
        check_shared_references(val, name)
    end

    # Validate optional bounds.
    !isnothing(ub_bounds) && isnothing(ub) && throw(ArgumentError("ub_bounds provided but ub is nothing"))
    !isnothing(uq_bounds) && isnothing(uq) && throw(ArgumentError("uq_bounds provided but uq is nothing"))
    !isnothing(m_bounds) && isnothing(m) && throw(ArgumentError("m_bounds provided but m is nothing"))
    isnothing(ub) || check_bounds_keys(ub_bounds, ub, "ub_bounds")
    isnothing(uq) || check_bounds_keys(uq_bounds, uq, "uq_bounds")
    isnothing(m) || check_bounds_keys(m_bounds, m, "m_bounds")

    uq_base_keys = isnothing(uq) ? () : _uq_base_keys(uq)
    uq_links_local = isnothing(uq) ? NamedTuple() : _uq_links(uq)

    # Validate links
    if !isempty(uq_links_local)
        for (child, link) in pairs(uq_links_local)
            parent = link.dependson
            parent_val = get(uq, parent, nothing)
            is_source_link(parent_val) && error("Dependent source $child cannot depend on another dependent ($parent)")
            isnothing(parent_val) && error("Dependent source $child points to missing parent $parent in uq")
            child == parent && error("Dependent source $child cannot depend on itself")
        end
    end

    if !isnothing(Qₛ)
        for (k, v) in pairs(uq)
            if is_source_link(v) && !isnothing(Qₛ[k])
                throw(ArgumentError("Qₛ entry for dependent source $k must be nothing"))
            end
        end
    end
    if !isnothing(uq_bounds)
        for (k, v) in pairs(uq)
            if is_source_link(v) && !isnothing(uq_bounds[k])
                throw(ArgumentError("uq_bounds entry for dependent source $k must be nothing"))
            end
        end
    end

    # Optionally cache precomputed mass-fraction steps and water-mass matrix.
    m_steps = nothing
    if cache_mass_fraction_steps
        isnothing(m) && throw(ArgumentError("cache_mass_fraction_steps=true requires m"))
        m_steps = precompute_mass_fraction_steps(m, γ)
    end

    A_cached = nothing
    if cache_A
        isnothing(m) && throw(ArgumentError("cache_A=true requires m"))
        if cache_mass_fraction_steps
            A_cached = watermassmatrix(m, γ, m_steps)
        else
            A_cached = watermassmatrix(m, γ)
        end
    end

    vector = cache_vector ? vectorize_controls(ub, uq, m) : nothing
    ub_length = isnothing(ub) ? nothing : length(vec(ub))
    uq_length = isnothing(uq) ? nothing : length(_vec_base_sources(uq, uq_base_keys))
    m_length = isnothing(m) ? nothing : length(vec(m))

    dub_cached = isnothing(ub) ? nothing : deepcopy(ub)
    gdub_cached = isnothing(ub) ? nothing : deepcopy(ub)
    b_cached = isnothing(u₀) ? nothing : deepcopy(u₀)
    boundary = BoundaryControls(ub, u₀, Qᵤ, ub_bounds, dub_cached, gdub_cached, b_cached)

    duq_cached = isnothing(uq) ? nothing : _make_source_buffer(uq, q₀)
    gduq_cached = isnothing(uq) ? nothing : _make_source_buffer(uq, q₀)
    q_cached = isnothing(uq) ? (isnothing(q₀) ? nothing : deepcopy(q₀)) : _make_source_buffer(uq, q₀)
    source = SourceControls(uq, q₀, Qₛ, uq_bounds, duq_cached, gduq_cached, q_cached)

    gm_cached = isnothing(m) ? nothing : deepcopy(m)
    massfrac = MassFractionControls(m, m₀, Qₘ, m_bounds, m_steps, A_cached, gm_cached)

    return ControlParameters{typeof(boundary), typeof(source), typeof(massfrac),
                             typeof(vector), typeof(ub_length), typeof(uq_length),
                             typeof(m_length), typeof(γ)}(
        boundary, source, massfrac, vector, ub_length, uq_length, m_length, γ)
end

function vectorize_controls(ub::NamedTuple, uq::NamedTuple, m::NamedTuple)
    uq_base_keys = _uq_base_keys(uq)
    return vcat(vec(ub), _vec_base_sources(uq, uq_base_keys), vec(m))
end

function vec(cp::ControlParameters)
    if !isnothing(cp.vector)
        return cp.vector
    else
        return vectorize_controls(cp.boundary.ub, cp.source.uq, cp.massfrac.m)
    end
end

function unvec(cp::ControlParameters, vector::AbstractVector)
    u_len = cp.ub_length
    q_len = cp.uq_length
    m_len = cp.m_length
    idx = 1
    ub = nothing; uq = nothing; m = nothing
    if !isnothing(cp.boundary.ub)
        ub = unvec(cp.boundary.ub, vector[idx:idx+u_len-1])
        idx += u_len
    end
    if !isnothing(cp.source.uq)
        uq = deepcopy(cp.source.uq)
        base_keys = _uq_base_keys(cp.source.uq)
        _unvec_sources!(uq, vector; idx = idx, base_keys = base_keys)
        idx += q_len
    end
    if !isnothing(cp.massfrac.m)
        m = unvec(cp.massfrac.m, vector[idx:idx+m_len-1])
    end
    return ub, uq, m
end

function unvec!(cp::ControlParameters, vector::AbstractVector)
    idx = 1
    if !isnothing(cp.boundary.ub)
        idx = unvec!(cp.boundary.ub, vector; idx = idx, return_idx = true)
    end
    if !isnothing(cp.source.uq)
        base_keys = _uq_base_keys(cp.source.uq)
        idx = _unvec_sources!(cp.source.uq, vector; idx = idx, return_idx = true, base_keys = base_keys)
    end
    if !isnothing(cp.massfrac.m)
        unvec!(cp.massfrac.m, vector; idx = idx)
    end
    return cp.boundary.ub, cp.source.uq, cp.massfrac.m
end

function _vectorize_bounds(bnds::Tuple{T,T}, controls::NamedTuple) where {T<:Real}
    lb = Float64[]
    ub = Float64[]
    lower, upper = bnds
    for control in controls
        if isnothing(control) || is_source_link(control)
            continue
        end
        n = length(vec(control))
        append!(lb, fill(T(lower), n))
        append!(ub, fill(T(upper), n))
    end
    return lb, ub
end

function _vectorize_bounds(bnds::NamedTuple, controls::NamedTuple)
    keys(bnds) == keys(controls) || error("bounds keys must match controls")
    lb = Float64[]
    ub = Float64[]
    for key in eachindex(controls)
        control = controls[key]
        isnothing(control) && continue
        is_source_link(control) && continue
        bnd = bnds[key]
        n = length(vec(control))
        if control isa NamedTuple
            if bnd isa Tuple
                lb_sub, ub_sub = _vectorize_bounds(bnd, control)
                append!(lb, lb_sub); append!(ub, ub_sub)
            elseif bnd isa NamedTuple
                keys(bnd) == keys(control) || error("nested bounds keys for $key must match control vector keys")
                lb_sub, ub_sub = _vectorize_bounds(bnd, control)
                append!(lb, lb_sub); append!(ub, ub_sub)
            end
        else
            lower, upper = bnd
            (lower isa Real && upper isa Real) || error("bounds for $key must be scalars")
            append!(lb, fill(Float64(lower), n)); append!(ub, fill(Float64(upper), n))
        end
    end
    return lb, ub
end

function vectorize_control_parameter_bounds(cp::ControlParameters)
    total_len = 0
    total_len += isnothing(cp.ub_length) ? 0 : cp.ub_length
    total_len += isnothing(cp.uq_length) ? 0 : cp.uq_length
    total_len += isnothing(cp.m_length) ? 0 : cp.m_length
    lower_bounds = Vector{Float64}(undef, total_len)
    upper_bounds = Vector{Float64}(undef, total_len)
    idx = 1
    if !isnothing(cp.boundary.ub)
        lower_ub, upper_ub = _vectorize_bounds(cp.boundary.ub_bounds, cp.boundary.ub)
        len = length(lower_ub)
        lower_bounds[idx:idx+len-1] .= lower_ub
        upper_bounds[idx:idx+len-1] .= upper_ub
        idx += len
    end
    if !isnothing(cp.source.uq)
        lower_uq, upper_uq = _vectorize_bounds(cp.source.uq_bounds, cp.source.uq)
        len = length(lower_uq)
        lower_bounds[idx:idx+len-1] .= lower_uq
        upper_bounds[idx:idx+len-1] .= upper_uq
        idx += len
    end
    if !isnothing(cp.massfrac.m)
        lower_m, upper_m = _vectorize_bounds(cp.massfrac.m_bounds, cp.massfrac.m)
        len = length(lower_m)
        lower_bounds[idx:idx+len-1] .= lower_m
        upper_bounds[idx:idx+len-1] .= upper_m
    end
    return lower_bounds, upper_bounds
end

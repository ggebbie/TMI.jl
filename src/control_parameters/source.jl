is_source_link(x) = x isa NamedTuple && haskey(x, :scale) && haskey(x, :dependson)

function _uq_base_keys(uq::NamedTuple)
    bases = Symbol[]
    for (k, v) in pairs(uq)
        if !is_source_link(v) && !isnothing(v)
            push!(bases, k)
        end
    end
    return Tuple(bases)
end

function _uq_links(uq::NamedTuple)
    links = Pair{Symbol, Any}[]
    for (k, v) in pairs(uq)
        if is_source_link(v)
            push!(links, k => (dependson = Symbol(v.dependson), scale = v.scale))
        end
    end
    return isempty(links) ? NamedTuple() : (; links...)
end

function _make_source_buffer(uq::NamedTuple, q₀::Union{Nothing, NamedTuple})
    keys_uq = keys(uq)
    vals = Vector{Any}(undef, length(uq))
    for (i, key) in enumerate(keys_uq)
        v = uq[key]
        if is_source_link(v)
            q0_val = isnothing(q₀) ? nothing : get(q₀, key, nothing)
            if isnothing(q0_val)
                parent = Symbol(v.dependson)
                parent_val = get(uq, parent, nothing)
                isnothing(parent_val) && error("Dependent source $key requires parent $parent in uq")
                is_source_link(parent_val) && error("Dependent source $key cannot depend on another dependent ($parent)")
                q0_val = deepcopy(parent_val)
                zero!(q0_val)
            end
            vals[i] = deepcopy(q0_val)
        else
            vals[i] = deepcopy(v)
        end
    end
    return NamedTuple{keys_uq}(vals)
end

function _vec_base_sources(uq::NamedTuple, base_keys::Tuple)
    res = Float64[]
    for key in base_keys
        v = uq[key]
        isnothing(v) && continue
        append!(res, vec(v))
    end
    return res
end

function _unvec_sources!(uq::NamedTuple, uvec::AbstractVector;
                         idx::Int = 1,
                         return_idx::Bool = false,
                         base_keys::Tuple = _uq_base_keys(uq))
    for key in base_keys
        v = uq[key]
        if !isnothing(v)
            idx = unvec!(v, uvec; idx = idx, return_idx = true)
        end
    end
    return return_idx ? idx : nothing
end

struct SourceControls{UQ,Q0,QS,UQB,DUQ,GDUQ,BQ}
    uq::UQ
    q₀::Q0
    Qₛ::QS
    uq_bounds::UQB
    duq::DUQ
    gduq::GDUQ
    q::BQ
end

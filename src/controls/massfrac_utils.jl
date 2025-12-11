wet(m::MassFraction) = m.γ.wet

function Base.similar(v::Vector{<:MassFraction})
    # return MassFraction[similar(mf) for mf in v]
    return map(similar, v)
end

function Base.similar(nt::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    return map(similar, nt)
end

function Base.similar(mf::MassFraction)
    similar_m_fraction = similar(mf.fraction)
    similar_m_fraction[.!(mf.γ.interior)] .= NaN

    return MassFraction(
        similar_m_fraction,
        mf.γ,
        mf.name,
        mf.longname,
        mf.units,
        mf.position
    )
end

function zero!(m::MassFraction)
    m.fraction[m.γ.wet] .= 0
    return m
end

function zero!(m::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    for v in m
        zero!(v)
    end
    return m
end

function Base.vec(u::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}) where names
    if isempty(u)
        return Float64[]
    end
    T = eltype(values(u)[1].fraction)
    #T = eltype(u)
    uvec = Vector{T}(undef,0)
    for v in u
        #append!(uvec,v.tracer[v.wet])
        append!(uvec,vec(v))
    end
    return uvec
end

function unvec!(u::MassFraction, uvec::Vector{T}; idx::Int = 1, return_idx::Bool = false) where T <: Real
    mask = wet(u)
    data = u.fraction
    @inbounds for I in eachindex(mask)
        if mask[I]
            data[I] = uvec[idx]
            idx += 1
        end
    end
    return return_idx ? idx : nothing
end

function unvec!(u::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}, uvec::Vector; idx::Int = 1, return_idx::Bool = false) where names
    for v in u
        idx = unvec!(v, uvec; idx = idx, return_idx = true)
    end
    return return_idx ? idx : nothing
end

function adjustmassfraction!(m::MassFraction{T}, uvec::AbstractVector{T};
                             idx::Int = 1,
                             return_idx::Bool = false,
                             r::Real = 1.0) where {T<:Real}
    mask = wet(m)
    data = m.fraction
    @inbounds for I in eachindex(mask)
        if mask[I]
            data[I] += r .* uvec[idx]
            idx += 1
        end
    end
    return return_idx ? idx : nothing
end

function adjustmassfraction!(m::NamedTuple{names, <:Tuple{Vararg{MassFraction}}}, uvec::AbstractVector;
                             idx::Int = 1,
                             return_idx::Bool = false,
                             r::Real = 1.0) where {names}
    for v in m
        idx = adjustmassfraction!(v, uvec; idx = idx, return_idx = true, r = r)
    end
    return return_idx ? idx : nothing
end

function unvec(u₀::NamedTuple{names, <:Tuple{Vararg{MassFraction}}},uvec::Vector) where names
    u = deepcopy(u₀)
    unvec!(u,uvec)
    return u
end


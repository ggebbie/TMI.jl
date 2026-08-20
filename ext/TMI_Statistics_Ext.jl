module TMI_Statistics_Ext

using TMI, Statistics
import TMI: zonalaverage

"""
    mean(field::Field)
    mean(field::Field; dims)
    mean(boundary::BoundaryCondition, γ::Grid; dims)
    zonalaverage(field::Field)
    zonalaverage(boundary::BoundaryCondition, γ::Grid)

Take a physically weighted mean, optionally removing longitude.

# Arguments
- `field`: gridded quantity
- `boundary`, `γ`: surface quantity and its grid
- `dims`: dimensions to average

# Output
- `average`: weighted scalar or array
"""
function Statistics.mean(field::Field; dims=nothing)
    volume = cellvolume(field.γ)
    if isnothing(dims)
        weights = volume.tracer[wet(field)]
        return sum(field.tracer[wet(field)] .* weights) / sum(weights)
    end
    valid = wet(field)
    values, weights = copy(field.tracer), copy(volume.tracer)
    values[.!valid] .= zero(eltype(values))
    weights[.!valid] .= zero(eltype(weights))
    return sum(values .* weights; dims) ./ sum(weights; dims)
end

function Statistics.mean(boundary::BoundaryCondition, γ::Grid; dims)
    area = cellarea(γ)
    valid = boundary.wet
    values, weights = copy(boundary.tracer), copy(area.tracer)
    values[.!valid] .= zero(eltype(values))
    weights[.!valid] .= zero(eltype(weights))
    return sum(values .* weights; dims) ./ sum(weights; dims)
end

function zonalaverage(field::Field)
    longitude = findfirst(axis -> axis === field.γ.lon, field.γ.axes)
    return dropdims(mean(field; dims=longitude); dims=longitude)
end

function zonalaverage(boundary::BoundaryCondition, γ::Grid)
    longitude = findfirst(axis -> axis == γ.lon, boundary.axes)
    return dropdims(mean(boundary, γ; dims=longitude); dims=longitude)
end

end

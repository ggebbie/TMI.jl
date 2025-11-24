"""
    Observations{T, V, W}

Container for observation data with precomputed interpolation weights.
"""
struct Observations{T, V, W}
    values::Vector{T}
    locs::Vector{V}
    wis::Vector{W}
end

"""
    Observations(c_obs, locs, γ) -> Observations

Construct Observations from values, locations, and grid; precomputes interpolation weights.
"""
function Observations(values::Vector{T}, locs::Vector{V}, γ::Grid) where {T, V}
    length(values) == length(locs) || error("values and locs must have the same length")
    wis = [interpindex(loc, γ) for loc in locs]
    return Observations(values, locs, wis)
end

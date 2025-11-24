"""
    Observations{T, V, W}

Container for observation data with precomputed interpolation weights.
"""
struct Observations{T, V, W}
    c_obs::Vector{T}
    locs::Vector{V}
    wis::Vector{W}
end

"""
    Observations(c_obs, locs, γ) -> Observations

Construct Observations from values, locations, and grid; precomputes interpolation weights.
"""
function Observations(c_obs::Vector{T}, locs::Vector{V}, γ::Grid) where {T, V}
    length(c_obs) == length(locs) || error("c_obs and locs must have the same length")
    wis = [interpindex(loc, γ) for loc in locs]
    return Observations(c_obs, locs, wis)
end

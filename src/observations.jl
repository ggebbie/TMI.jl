"""
    Observations{T, V, WI}

Container for observation data with precomputed interpolation weights and optional weighting matrix.
Stores observational data along with their spatial locations and interpolation information for efficient
comparison with model fields. Multiple Observations can be organized in a NamedTuple for multi-tracer problems.

# Type Parameters
- `T`: Type of observation values (typically Float64)
- `V`: Type of location data (e.g., Tuple of coordinates)
- `WI`: Type of interpolation weight indices

# Fields
- `values`: Vector of observed values
- `locs`: Vector of spatial locations corresponding to each observation
- `wis`: Vector of precomputed interpolation indices/weights for efficient model-data comparison
- `W`: Optional inverse covariance matrix (Symmetric) for observation errors

# Constructor
```julia
Observations(values::Vector, locs::Vector, γ::Grid; W = nothing)
```

Constructs an Observations object from observation values, their locations, and a computational grid.
Automatically precomputes interpolation weights for efficient model-data interpolation.

# Arguments
- `values`: Vector of observation values
- `locs`: Vector of observation locations (must have same length as `values`)
- `γ`: Grid object defining the computational domain
- `W`: Optional inverse covariance matrix (Symmetric) for observation uncertainties

# Examples
```julia
# Basic usage without weighting
obs = Observations([1.0, 2.0, 3.0], [loc1, loc2, loc3], γ)

# With inverse covariance matrix for uncertainty weighting
W_matrix = Symmetric(inv(cov_matrix))
obs_θ = Observations(θ_values, θ_locs, γ; W = W_matrix)

# Multiple tracers organized in a NamedTuple
obs = (θ = obs_θ, S = obs_S)
```

# Notes
The interpolation weights (`wis`) are computed at construction time using `interpindex(loc, γ)` for
each location, enabling efficient repeated interpolation operations during optimization.
"""
struct Observations{T, V, WI}
    values::Vector{T}
    locs::Vector{V}
    wis::Vector{WI}
    W::Union{Symmetric, Nothing}
end

"""
    Observations(values, locs, γ; W = nothing) -> Observations

Construct Observations from values, locations, and grid; precomputes interpolation weights.
Optional weighting matrix W can be provided.
"""
function Observations(values::Vector{T}, locs::Vector{V}, γ::Grid; W::Union{Symmetric, Nothing} = nothing) where {T, V}
    length(values) == length(locs) || error("values and locs must have the same length")

    # Validate weighting matrix dimensions if provided
    if !isnothing(W)
        n = length(values)
        size(W) == (n, n) || error("W must have dimensions ($n, $n) to match the length of values")
    end

    wis = [interpindex(loc, γ) for loc in locs]
    return Observations{T, V, typeof(first(wis))}(values, locs, wis, W)
end

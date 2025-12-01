"""
    Observations{T, V, WI}
    Observations(values; locs=nothing, γ=nothing, W=nothing)

Container for observation data with optional locations and weighting. If locations are provided,
interpolation weights are precomputed via `interpindex(loc, γ)`.

- `values::Union{Vector, Field}`: observation values (vector or grid Field)
- `locs::Union{Vector, Nothing}`: observation locations (or `nothing`)
- `γ::Grid`: grid used for interpolation (required when `locs` are given)
- `W::Union{Symmetric, Diagonal, Nothing}`: optional weighting matrix
"""
struct Observations{T, V, WI}
    values::Union{Vector{T}, Field{T}}
    locs::Union{Vector{V}, Nothing}
    wis::Union{Vector{WI}, Nothing}
    W::Union{Symmetric, Diagonal, Nothing}
end

"""
    Observations(values; locs = nothing, γ = nothing, W = nothing) -> Observations

Construct Observations from values with optional locations, grid, and weighting matrix.
If locations are provided, interpolation weights are precomputed.
"""
function Observations(values::Union{Vector{T}, Field{T}};
                      locs::Union{Vector{V}, Nothing} = nothing,
                      γ::Union{Grid, Nothing} = nothing,
                      W::Union{Symmetric, Diagonal, Nothing} = nothing) where {T, V}
    if !isnothing(W)
        # Keep a quick sanity check that weights match the observation count.
        n = values isa Vector ? length(values) : length(vec(values))
        size(W) == (n, n) || error("W must have dimensions ($n, $n) to match the length of values")
    end

    if isnothing(locs)
        wis = nothing
        VType = Nothing
        WIType = Nothing
    else
        values isa Field && error("locs should be `nothing` when values is a Field")
        length(values) == length(locs) || error("values and locs must have the same length")
        isnothing(γ) && error("γ must be provided when locs are given")
        # Precompute interpolation weights once; reused by observe/gobserve.
        wis = [interpindex(loc, γ) for loc in locs]
        VType = V
        WIType = eltype(wis)
    end

    return Observations{T, VType, WIType}(values, locs, wis, W)
end

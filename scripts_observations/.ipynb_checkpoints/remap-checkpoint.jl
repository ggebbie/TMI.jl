"""
    remap_fine_to_coarse(
        lon_coarse, lat_coarse,
        wetmask_coarse::Matrix{Bool},
        lon_fine, lat_fine,
        var_fine::Matrix{<:Real};
        weights_fine::Union{Nothing, Matrix{<:Real}} = nothing,
        method::Union{Symbol, Function} = :mean
    ) -> Matrix{Union{Missing, Float64}}

Coarsens a fine-grid variable `var_fine` (on lon_fine × lat_fine) to a coarser
lon_coarse × lat_coarse grid using the provided `wetmask_coarse` to restrict output
to ocean (true) cells only.

NaNs in `var_fine` are ignored. Supports weighted or unweighted averaging or summing.

# Arguments
- `lon_coarse`, `lat_coarse`: 1D vectors defining the coarse grid centers
- `wetmask_coarse`: matrix of shape `(nlat_coarse, nlon_coarse)` with `true` for ocean, `false` for land
- `lon_fine`, `lat_fine`: 1D vectors defining the fine grid centers
- `var_fine`: 2D matrix of shape `(nlat_fine, nlon_fine)` containing the fine-grid variable values

# Keyword Arguments
- `weights_fine`: optional matrix of weights with shape `(nlat_fine, nlon_fine)`
- `method`: aggregation strategy (`:mean` [default], `:sum`, or a custom function)

# Returns
- A matrix of shape `(nlat_coarse, nlon_coarse)` with remapped values, or `missing` for land cells.

# Notes
All 2D matrices must follow Julia's default memory layout:
`(latitude index, longitude index)` → that is, **rows = lat, columns = lon**.
"""
function remap_fine_to_coarse(
    lon_coarse::AbstractVector, lat_coarse::AbstractVector,
    wetmask_coarse::AbstractMatrix{Bool}, γ_coarse::Grid, 
    lon_fine::AbstractVector, lat_fine::AbstractVector,
    var_fine::Matrix{<:Real};
    weights_fine::Union{Nothing, Matrix{<:Real}} = nothing,
    method::Union{Symbol, Function} = :mean
)
    nlat2, nlon2 = length(lat_coarse), length(lon_coarse)
    nlat1, nlon1 = length(lat_fine), length(lon_fine)

    @assert size(wetmask_coarse) == (nlat2, nlon2)
    @assert size(var_fine) == (nlat1, nlon1)
    if weights_fine !== nothing
        @assert size(weights_fine) == size(var_fine)
    end

    # Compute grid spacing (assumes uniform spacing)
    dlat2 = mean(diff(lat_coarse))
    dlon2 = mean(diff(lon_coarse))

    # Output: missing for land
    var_coarse = fill(NaN, nlat2, nlon2)  # type: Matrix{Float64}

    # Fine grid mesh
    lon1_grid = repeat(lon_fine', nlat1, 1)
    lat1_grid = repeat(lat_fine, 1, nlon1)

    lon_wrap(lon) = mod(lon, 360)

    for (i2, lat2) in enumerate(lat_coarse), (j2, lon2) in enumerate(lon_coarse)
        if !wetmask_coarse[i2, j2]
            continue
        end

        # Define bounds using half spacing
        lat_min, lat_max = lat2 - dlat2/2, lat2 + dlat2/2
        lon_min, lon_max = lon_wrap(lon2 - dlon2/2), lon_wrap(lon2 + dlon2/2)

        lon1_wrapped = lon_wrap.(lon1_grid)

        in_lat = (lat1_grid .>= lat_min) .& (lat1_grid .<= lat_max)
        in_lon = if lon_min <= lon_max
            (lon1_wrapped .>= lon_min) .& (lon1_wrapped .<= lon_max)
        else
            (lon1_wrapped .>= lon_min) .| (lon1_wrapped .<= lon_max)
        end

        region_mask = in_lat .& in_lon

        # Extract valid values
        vals_all = vec(var_fine[region_mask])
        valid_idx = .!isnan.(vals_all)
        vals = vals_all[valid_idx]

        if isempty(vals)
            continue
        end

        if weights_fine !== nothing
            w_all = vec(weights_fine[region_mask])
            w = w_all[valid_idx]

            if method == :mean
                var_coarse[i2, j2] = sum(vals .* w) / sum(w)
            elseif method == :sum
                var_coarse[i2, j2] = sum(vals .* w)
            elseif isa(method, Function)
                var_coarse[i2, j2] = method(vals, w)
            else
                error("Unsupported method with weights: $method")
            end
        else
            if method == :mean
                var_coarse[i2, j2] = mean(vals)
            elseif method == :sum
                var_coarse[i2, j2] = sum(vals)
            elseif isa(method, Function)
                var_coarse[i2, j2] = method(vals)
            else
                error("Unsupported method: $method")
            end
        end
    end
    wetmask_new = (!).(isnan).(var_coarse')
    
    return BoundaryCondition(
        var_coarse',
        (lon_coarse, lat_coarse),
        γ_coarse.depth[1],
        3,
        1,
        wetmask_new,
        :nothing,
        "coarsened ",
        "nothing"
    )
    
    # return var_coarse
end

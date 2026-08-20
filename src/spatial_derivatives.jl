"""
    adjacencymatrix(γ, mask, offsets)
    surfaceadjacencymatrix(γ, offsets)
    interioradjacencymatrix(γ, offsets)

Build sparse adjacency matrices for selected, surface, or interior cells.

# Arguments
- `γ::Grid`: grid and wrapping rules
- `mask`: cells included as rows and columns
- `offsets`: relative neighbor directions

# Output
- `G`: sparse adjacency matrix in the selected wet-cell order
"""
function adjacencymatrix(
    γ::Grid{R,N},
    mask::AbstractArray{Bool,N},
    offsets::NTuple{M,CartesianIndex{N}},
) where {R,N,M}
    cells = cartesianindex(mask)
    row_at = zeros(Int, size(mask))
    row_at[mask] = eachindex(cells)
    rows, columns = Int[], Int[]
    for (row, cell) in enumerate(cells)
        for offset in offsets
            neighbor, inbounds = step_cartesian(cell, offset, γ)
            inbounds && mask[neighbor] || continue
            push!(rows, row)
            push!(columns, row_at[neighbor])
        end
    end
    return sparse(rows, columns, ones(length(rows)), length(cells), length(cells))
end
function surfaceadjacencymatrix(
    γ::Grid{R,3},
    offsets::NTuple{M,CartesianIndex{3}},
) where {R,M}
    surface = surfaceindex(γ)
    mask = falses(size(γ.wet))
    mask[:, :, surface] .= γ.wet[:, :, surface]
    return adjacencymatrix(γ, mask, offsets)
end
function interioradjacencymatrix(
    γ::Grid{R,N},
    offsets::NTuple{M,CartesianIndex{N}},
) where {R,N,M}
    return adjacencymatrix(γ, γ.interior, offsets)
end

"""
    surfacelaplacianmatrix(γ)

Construct the horizontal surface Laplacian with wrapping and zero land flux.

# Arguments
- `γ::Grid`: three-dimensional TMI grid

# Output
- `L`: sparse wet-cell Laplacian with units km⁻²
"""
function surfacelaplacianmatrix(γ::Grid{R,3}) where R
    surface = surfaceindex(γ)
    cells = findall(γ.wet[:, :, surface])
    zonal_adjacency = surfaceadjacencymatrix(
        γ,
        (CartesianIndex(-1, 0, 0), CartesianIndex(1, 0, 0)),
    )
    meridional_adjacency = surfaceadjacencymatrix(
        γ,
        (CartesianIndex(0, -1, 0), CartesianIndex(0, 1, 0)),
    )
    dx = zonalgriddist(γ) ./ 1000
    dy = haversine((γ.lon[1], γ.lat[1]), (γ.lon[1], γ.lat[2])) / 1000
    inverse_dx² = spdiagm(0 => [inv(dx[cell[2]]^2) for cell in cells])
    inverse_dy² = spdiagm(0 => fill(inv(dy^2), length(cells)))
    weighted_adjacency = inverse_dx² * zonal_adjacency +
        inverse_dy² * meridional_adjacency
    return weighted_adjacency -
        spdiagm(0 => vec(sum(weighted_adjacency; dims=2)))
end

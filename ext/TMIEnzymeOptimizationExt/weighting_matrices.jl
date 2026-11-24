# Weighting matrix of the prior smoothness cost of boundary adjustments: the
# Laplacian on the wet surface cells (`surfacelaplacianmatrix`).

"""
    function adjacencymatrix(γ, mask, offsets)

Sparse matrix with a 1 where the cell of column j is a neighbor of the cell of
row i. Rows and columns follow the order of the cells in `mask`.

# Arguments
- `γ::Grid`: TMI grid, which sets zonal wrapping
- `mask`: cells included as rows and columns
- `offsets`: steps to the neighbors, e.g. `(CartesianIndex(1, 0, 0),)`
# Output
- `G`: sparse adjacency matrix
"""
function adjacencymatrix(γ::Grid{T,N}, mask::BitArray{N},
    offsets::NTuple{M,CartesianIndex{N}}) where {T,N,M}
    cells = cartesianindex(mask)
    R = linearindex(mask)
    rows, columns = Int[], Int[]
    for (row, cell) in enumerate(cells), offset in offsets
        neighbor, inbounds = step_cartesian(cell, offset, γ)
        inbounds && mask[neighbor] || continue
        push!(rows, row)
        push!(columns, R[neighbor])
    end
    return sparse(rows, columns, ones(length(rows)), length(cells), length(cells))
end

"""
    function surfacelaplacianmatrix(γ)

Horizontal Laplacian on wet surface cells, with zonal wrapping and no flux
across land. Rows and columns follow the order of `vec` of a surface
`BoundaryCondition`.

# Arguments
- `γ::Grid`: TMI grid
# Output
- `∇²`: sparse Laplacian [km⁻²]
"""
function surfacelaplacianmatrix(γ::Grid{R,3}) where R
    surface = surfaceindex(γ)
    mask = falses(size(γ.wet))
    mask[:, :, surface] .= γ.wet[:, :, surface]
    zonal = adjacencymatrix(γ, mask, (CartesianIndex(-1, 0, 0), CartesianIndex(1, 0, 0)))
    meridional = adjacencymatrix(γ, mask, (CartesianIndex(0, -1, 0), CartesianIndex(0, 1, 0)))
    dx = TMI.zonalgriddist(γ) ./ 1000
    dy = TMI.haversine((γ.lon[1], γ.lat[1]), (γ.lon[1], γ.lat[2])) / 1000
    cells = cartesianindex(mask)
    weighted = spdiagm(0 => [inv(dx[cell[2]]^2) for cell in cells]) * zonal +
        meridional / dy^2
    return weighted - spdiagm(0 => vec(sum(weighted; dims=2)))
end

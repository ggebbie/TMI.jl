"""
    euclidean_distance(latᵢ, lonᵢ, depthᵢ, latⱼ, lonⱼ, depthⱼ; R=6371.0)

Straight-line distance between two points below the sea surface, accounting
for their depths below a sphere of radius `R`.

# Arguments
- `latᵢ`, `lonᵢ`, `depthᵢ`: latitude [°N], longitude [°E], and depth [m] of the first point
- `latⱼ`, `lonⱼ`, `depthⱼ`: latitude [°N], longitude [°E], and depth [m] of the second point
- `R`: radius of the sphere [km]

# Output
- `d`: distance [km]
"""
function euclidean_distance(latᵢ, lonᵢ, depthᵢ, latⱼ, lonⱼ, depthⱼ; R=6_371.0)
    ϕᵢ, λᵢ, ϕⱼ, λⱼ = deg2rad.((latᵢ, lonᵢ, latⱼ, lonⱼ))
    rᵢ, rⱼ = R - depthᵢ / 1_000, R - depthⱼ / 1_000
    xyz(r, ϕ, λ) = (r * cos(ϕ) * cos(λ), r * cos(ϕ) * sin(λ), r * sin(ϕ))
    return sqrt(sum(abs2, xyz(rⱼ, ϕⱼ, λⱼ) .- xyz(rᵢ, ϕᵢ, λᵢ)))
end

"""
    adjacencymatrix(γ, mask, offsets)

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
    gaussianprecision(σ, L, γ)

Sparse inverse covariance of gridded interior observations with a Gaussian
horizontal correlation, built by a Vecchia approximation on horizontal neighbors.

# Arguments
- `σ`: standard deviation at each interior cell
- `L`: decorrelation length scale [km] at each interior cell
- `γ`: TMI grid

# Output
- `Wⁱ`: symmetric sparse inverse covariance matrix
"""
function gaussianprecision(σ::AbstractVector, L::AbstractVector, γ::Grid{R,3}) where R
    T = float(promote_type(eltype(σ), eltype(L), R))
    σ, L = T.(σ), T.(L)
    cells = cartesianindex(γ.interior)
    horizontal = (CartesianIndex(-1, 0, 0), CartesianIndex(1, 0, 0),
        CartesianIndex(0, -1, 0), CartesianIndex(0, 1, 0))
    rows, columns, _ = findnz(tril(adjacencymatrix(γ, γ.interior, horizontal), -1))
    neighbors = [Int[] for _ in cells]
    foreach((i, j) -> push!(neighbors[i], j), rows, columns)

    function correlation(i, j)
        i == j && return one(T)
        ci, cj = cells[i], cells[j]
        d = euclidean_distance(γ.lat[ci[2]], γ.lon[ci[1]], γ.depth[ci[3]],
            γ.lat[cj[2]], γ.lon[cj[1]], γ.depth[cj[3]])
        L² = L[i]^2 + L[j]^2
        return (2L[i] * L[j] / L²)^(3 / 2) * exp(-2d^2 / L²)
    end

    B_rows, B_columns, B_values = Int[], Int[], T[]
    conditional_variances = ones(T, length(cells))
    for i in eachindex(cells)
        push!(B_rows, i); push!(B_columns, i); push!(B_values, one(T))
        isempty(neighbors[i]) && continue
        RNN = [correlation(j, k) for j in neighbors[i], k in neighbors[i]]
        rNi = [correlation(j, i) for j in neighbors[i]]
        coefficients = cholesky(Hermitian(RNN)) \ rNi
        conditional_variances[i] = 1 - dot(rNi, coefficients)
        append!(B_rows, fill(i, length(coefficients)))
        append!(B_columns, neighbors[i])
        append!(B_values, -coefficients)
    end
    B = sparse(B_rows, B_columns, B_values, length(cells), length(cells))
    Dσ⁻¹ = spdiagm(0 => inv.(σ))
    Wⁱ = Dσ⁻¹ * transpose(B) * spdiagm(0 => inv.(conditional_variances)) * B * Dσ⁻¹
    return (Wⁱ + transpose(Wⁱ)) / 2
end

"""
    surfacelaplacianmatrix(γ)

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

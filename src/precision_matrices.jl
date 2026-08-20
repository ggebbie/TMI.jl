"""
    euclidean_distance(latᵢ, lonᵢ, depthᵢ, latⱼ, lonⱼ, depthⱼ; R)

Return the three-dimensional Euclidean distance in km between two geographic
locations, accounting for their depths below a spherical surface of radius `R`.

# Arguments
- `latᵢ`, `lonᵢ`, `depthᵢ`, `latⱼ`, `lonⱼ`, `depthⱼ`: geographic coordinates
- `R`: spherical radius in km

# Output
- `distance`: straight-line distance in km
"""
function euclidean_distance(
    latᵢ::Real, lonᵢ::Real, depthᵢ::Real,
    latⱼ::Real, lonⱼ::Real, depthⱼ::Real;
    R::Real=6_371.0,
)
    ϕᵢ, λᵢ, ϕⱼ, λⱼ = deg2rad.((latᵢ, lonᵢ, latⱼ, lonⱼ))
    rᵢ, rⱼ = R - depthᵢ / 1_000, R - depthⱼ / 1_000
    xyz(r, ϕ, λ) = (r * cos(ϕ) * cos(λ), r * cos(ϕ) * sin(λ), r * sin(ϕ))
    xᵢ, xⱼ = xyz(rᵢ, ϕᵢ, λᵢ), xyz(rⱼ, ϕⱼ, λⱼ)
    return sqrt(sum(abs2, xⱼ .- xᵢ))
end

"""
    gaussianprecision(σ, L)
    gaussianprecision(σ, L, γ)

Construct a sparse horizontal Vecchia precision for gridded observations.

# Arguments
- `σ`: positive pointwise standard deviations
- `L`: positive horizontal decorrelation scales in km
- `γ::Grid`: grid defining vector order

# Output
- `Q`: symmetric sparse inverse covariance
"""
function gaussianprecision(σ::Field{T,R,3}, L::Field) where {T,R}
    return gaussianprecision(σ.tracer[σ.γ.interior],
        L.tracer[σ.γ.interior], σ.γ)
end
function gaussianprecision(σ::AbstractVector{T}, L::AbstractVector{S},
    γ::Grid{R,3}) where {T<:Real,S<:Real,R<:Real}
    calculation_type = float(promote_type(T, S, R))
    σ_vector = calculation_type.(σ)
    L_vector = calculation_type.(L)

    cells = cartesianindex(γ.interior)
    adjacency = tril(interioradjacencymatrix(γ, (
        CartesianIndex(-1, 0, 0), CartesianIndex(1, 0, 0),
        CartesianIndex(0, -1, 0), CartesianIndex(0, 1, 0))), -1)
    neighbors = [Int[] for _ in cells]
    rows, columns, _ = findnz(adjacency)
    for (i, j) in zip(rows, columns)
        push!(neighbors[i], j)
    end

    function correlation(i, j)
        i == j && return one(calculation_type)
        ci, cj = cells[i], cells[j]
        d = euclidean_distance(γ.lat[ci[2]], γ.lon[ci[1]], γ.depth[ci[3]],
            γ.lat[cj[2]], γ.lon[cj[1]], γ.depth[cj[3]])
        Li, Lj = L_vector[i], L_vector[j]
        L² = Li^2 + Lj^2
        return (2Li * Lj / L²)^(3 / 2) * exp(-2d^2 / L²)
    end

    B_rows, B_columns, B_values = Int[], Int[], calculation_type[]
    conditional_variances = similar(L_vector)
    for i in eachindex(cells)
        push!(B_rows, i); push!(B_columns, i); push!(B_values, one(calculation_type))
        isempty(neighbors[i]) && (conditional_variances[i] = one(calculation_type); continue)
        RNN = [correlation(j, k) for j in neighbors[i], k in neighbors[i]]
        rNi = [correlation(j, i) for j in neighbors[i]]
        coefficients = cholesky(Hermitian(RNN)) \ rNi
        conditional_variances[i] = 1 - dot(rNi, coefficients)
        for (j, coefficient) in zip(neighbors[i], coefficients)
            push!(B_rows, i); push!(B_columns, j); push!(B_values, -coefficient)
        end
    end
    B = sparse(B_rows, B_columns, B_values, length(cells), length(cells))
    Dσ⁻¹ = spdiagm(0 => inv.(σ_vector))
    D⁻¹ = spdiagm(0 => inv.(conditional_variances))
    Q = Dσ⁻¹ * transpose(B) * D⁻¹ * B * Dσ⁻¹
    return (Q + transpose(Q)) / 2
end

"""
    boundary_smoothness_precision_matrix(σ, L, laplacian)

Construct a symmetric precision penalizing boundary-anomaly roughness.

# Arguments
- `σ`: boundary standard deviations
- `L`: boundary decorrelation scales in km
- `laplacian`: sparse surface Laplacian

# Output
- `Q`: symmetric sparse boundary-smoothness precision
"""
function boundary_smoothness_precision_matrix(
    σ::BoundaryCondition,
    L::BoundaryCondition,
    laplacian::SparseMatrixCSC,
)
    weights = vec(L) .^ 4 ./ vec(σ) .^ 2
    Q = transpose(laplacian) * spdiagm(0 => weights) * laplacian
    return (Q + transpose(Q)) / 2
end

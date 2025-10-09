"""
    block_sum(A, di, dj)

Collapse `A` into non–overlapping `di×dj` blocks by summing each block,
ignoring any `NaN` values within the block.

# Arguments
- `A::AbstractMatrix{T}` – fine‐grid array (dimensions must be divisible by `di` and `dj`)
- `di::Int` – block size in the row direction
- `dj::Int` – block size in the column direction
"""
function block_sum(A::AbstractMatrix{T}, di::Int, dj::Int) where T <: AbstractFloat
    m, n   = size(A)
    ni, nj = div(m, di), div(n, dj)
    B      = zeros(T, ni, nj)
    for i in 1:ni, j in 1:nj
        rows = (i-1)*di + 1 : i*di
        cols = (j-1)*dj + 1 : j*dj
        B[i, j] = sum(x for x in A[rows, cols] if !isnan(x))
    end
    return B
end

"""
    coarse_grain_boundary_condition(
        bc_fine::BoundaryCondition,
        γ_fine,
        γ_coarse,
        di::Int, dj::Int;
        method::Symbol = :average,
        use_area_weights::Bool = true,
        weights::Union{Nothing,AbstractMatrix{<:Real}} = nothing
    )

Coarse‐grain the fine‐grid BC `bc_fine` (on `γ_fine`) onto `γ_coarse` by
collapsing non‐overlapping `di×dj` blocks.

# Positional Arguments
- `bc_fine::BoundaryCondition`
- `γ_fine`       – fine grid (for computing cell areas)
- `γ_coarse`     – coarse grid (defines output coords & shape)
- `di, dj`       – block sizes in row/col directions

# Keyword Arguments
- `method::Symbol`            – `:sum` to return block‐sums, `:average` to return block‐averages (default)
- `use_area_weights::Bool`    – when `method==:average`, if `true` use area‐weighted average; if `false`, simple mean (default `true`)
- `weights`                   – if non‐`nothing`, use this fine‐grid weight array (masked by wet cells) instead of area or unit weights.

# Returns
A `BoundaryCondition` on `γ_coarse` with either the summed or averaged tracer.
"""
function coarse_grain_boundary_condition(
    bc_fine::BoundaryCondition,
    γ_fine::Grid,
    γ_coarse::Grid,
    di::Int, dj::Int;
    T::Type{<:AbstractFloat} = Float64,
    method::Symbol           = :average,
    use_area_weights::Bool   = true,
    weights::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
    var::Symbol = :nothing
)
    # 1) fine‐grid mask & area cast to T
    wet_fine  = T.(bc_fine.wet)                    # 1.0 where wet, 0.0 where dry
    area_fine = T.(surfacecellarea(γ_fine).tracer)        # always computed

    # 2) select fine‐grid weights W_fine (all T)
    W_fine = if weights !== nothing
        T.(weights)                     # user weights
    elseif use_area_weights
        area_fine .* wet_fine                             # area‐weighted
    else
        wet_fine                                         # simple count mask
    end
    Ttracer = T.(bc_fine.tracer)
    if var == :θ
        Ttracer .+= 273.15
    end
    W_fine[iszero.(W_fine)] .= NaN
    Ttracer[iszero.(Ttracer)] .= NaN
    # 3) block‐sum numerator & denominator
    num_coarse = block_sum(Ttracer .* W_fine, di, dj)
    den_coarse = block_sum(W_fine, di, dj)

    # 4) compute tracer_coarse per method
    tracer_coarse = similar(num_coarse, T)                # ensure element‐type T
    if method === :sum
        tracer_coarse .= num_coarse
    elseif method === :average
        @inbounds for idx in eachindex(num_coarse)
            tracer_coarse[idx] = den_coarse[idx] > zero(T) ?
                num_coarse[idx] / den_coarse[idx] :
                T(NaN)
        end
    else
        throw(ArgumentError("`method` must be :sum or :average"))
    end

    # 5) any positive weight ⇒ wet, else NaN in tracer_coarse
    wet_coarse = block_sum(wet_fine, di, dj) .> zero(T)
    tracer_coarse[.!wet_coarse] .= T(NaN)

    if var == :θ
        tracer_coarse .-= 273.15
    end
    # 6) construct new BoundaryCondition
    return BoundaryCondition(
        tracer_coarse,
        (γ_coarse.lon, γ_coarse.lat),
        γ_coarse.depth[1],
        3,
        1,
        wet_coarse,
        bc_fine.name,
        "coarsened " * bc_fine.longname,
        bc_fine.units
    )
end
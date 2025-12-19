
"""
    mass_conservation_constraints!(x, g; m0)

Fill the constraint vector `g` so each entry is the sum of the mass-fraction
controls at one interior grid cell. Pair this with equality bounds of `1.0` in
your solver to enforce that the mass fractions at every interior location sum
to one.

# Arguments
- `x`: Flat control vector. Contains updated mass fractions. 
- `g`: Preallocated constraint vector, length(g) == length(x).
- `m0`: Mass-fraction NamedTuple template used to reshape the relevant portion of `x`.

# Returns
- `g`: A constraint vector 
       populated in-place with the mass-fraction sums for each interior cell.
       Thus, should be a vector filled with ones
"""
function mass_conservation_constraints!(x::Vector{G}, g::Vector{G}; m0::NamedTuple) where G
    nm = length(vec(m0))
    m_x = x[end-nm+1:end]
    m = unvec(m0, m_x)

    firstm = first(m)
    γ = firstm.γ
    T = eltype(firstm.fraction)
    nint = sum(γ.interior)
    ngrid = size(γ.wet)

    total = zeros(T, ngrid)
    Iint = cartesianindex(γ.interior)
    for m1 in m
        sum(m1.γ.interior) == nint || error("mass fraction interior grids do not match")
    end
    for I in Iint
        s = zero(T)
        for m1 in m
            m1.γ.wet[I] || continue
            s += m1.fraction[I]
        end
        total[I] = s
    end

    total = total[firstm.γ.interior]    # interior ordering matches g
    for i in eachindex(g)               # solver-specific ordering for g
        g[i] = total[i]
    end
end
"""
    mass_conservation_jacobian!(x, rows, cols, values; m0)

Populate the Jacobian structure or values for
`mass_conservation_constraints!`. When `values === nothing`, the row/column
pattern is written. Otherwise the derivative entries are all ones in the entries 
corresponding to a mass fraction. 

# Arguments
- `x`: Control vector (only length is used to locate the mass-fraction block of the jacobian).
- `rows`, `cols`: Preallocated index arrays for the Jacobian in sparse format.
- `values`: Preallocated value array; pass `nothing` to fill only structure.
- `m0`: Mass-fraction template that defines the block layout and grid masks.

# Returns
- `rows`, `cols`, `values`: Filled in-place with the Jacobian sparsity and
  entries; returns `nothing`.
"""
function mass_conservation_jacobian!(x, rows, cols, values; m0)
    nm     = length(vec(m0))
    offset = length(x) - nm              # start of mass-fraction block in x

    firstm    = first(m0)
    interior  = firstm.γ.interior
    row_of_I  = zeros(Int, size(interior))
    r = 0
    for I in cartesianindex(interior)     # map each interior cell -> constraint row
        interior[I] || continue
        r += 1
        row_of_I[I] = r
    end

    col = offset
    nz  = 1
    oneT = one(eltype(x))
    for m1 in m0
        wet = m1.γ.wet
        for I in eachindex(wet)          # follow vec(m0) ordering (wet cells only)
            wet[I] || continue
            col += 1
            row = row_of_I[I]
            row == 0 && continue          # not an interior constraint -> skip
            if values === nothing
                rows[nz] = row
                cols[nz] = col
            else
                values[nz] = oneT         # ∂g_row/∂m_cell = 1
            end
            nz += 1
        end
    end
    return
end
"""
    foreachtracer(f, names)

The result of `f` for each tracer (or source) in `names`, as a `NamedTuple`
keyed by name.
"""
foreachtracer(f, names) = NamedTuple{names}(map(f, names))

"""
    Inversion(y; b₀, q₀=(;), r=(;), m₀, controls)

Inverse problem for adjustments `u` to first-guess boundary conditions `b₀`,
sources `q₀`, and mass fractions `m₀` that fit observations `y`, with cost
`J = nᵀ Wⁱ n + uvecᵀ Q⁻ uvec`.

# Arguments
- `y::Observations`: observations and their weighting matrix `Wⁱ`
- `b₀::NamedTuple`: first-guess surface `BoundaryCondition` for each tracer
- `q₀::NamedTuple`: first-guess interior `Source`s, linear scale
- `r::NamedTuple`: stoichiometric ratios, e.g. `(NO₃ = (qPO₄ = -15.5,),)`.
  Signs follow `A = watermassmatrix(m, γ)`, whose interior diagonal is +1
  (−1 in TMI files), so they are opposite to `r` in
  `steadyinversion(Alu, b, γ; q, r)` with a TMI file's `Alu`:
  `r = (PO₄ = (qPO₄ = -1.0,),)` gives phosphate from `qPO₄`.
- `m₀::NamedTuple`: first-guess `MassFraction`s
- `controls::NamedTuple`: adjusted quantities only, e.g.
  `(b = (θ = (σ=…, L=…, lb=…, ub=…),), q = (qPO₄ = (σ=…, lb=…, ub=…, logscale=true),), m = (σ=…, lb=…, ub=…))`.
  `σ` is the standard deviation of `u`, `L` [km] the smoothing scale of a
  boundary adjustment, and `lb`, `ub` are physical limits on `b`, `q`, and `m`;
  each is one number for every adjusted grid point or a field shaped like the
  first guess. A source with `logscale=true` is adjusted multiplicatively,
  `q = q₀ eᵘ`, so its `σ` is in log units (`σ = log(50)` is a one-σ factor of
  50), while its `lb` and `ub` remain limits on `q` itself.

# Output
- `Inversion`: first guesses, preallocated controls `u`, blocks `ranges` of
  `uvec`, and `σ`, `lb`, `ub`, `Q⁻` for `uvec`
"""
function Inversion(y::Observations; b₀::NamedTuple, q₀::NamedTuple=(;),
    r::NamedTuple=(;), m₀::NamedTuple, controls::NamedTuple)
    γ = y.γ
    bcontrols = get(controls, :b, (;))
    qcontrols = get(controls, :q, (;))
    mcontrols = get(controls, :m, nothing)

    # zero adjustments u for the controlled boundary conditions and sources; u.m
    # holds the adjusted mass fractions themselves, which watermassmatrix takes
    u = (b = foreachtracer(k -> zerosurfaceboundary(γ, b₀[k].name, b₀[k].longname, b₀[k].units),
            keys(bcontrols)),
        q = foreachtracer(k -> zerosource(γ, q₀[k].name, q₀[k].longname, q₀[k].units;
            logscale=get(qcontrols[k], :logscale, false)), keys(qcontrols)),
        m = isnothing(mcontrols) ? (;) : unvec(m₀, vec(m₀)))

    # controlled blocks of uvec, in b, q, m order
    blocks = [[(b₀[k], bcontrols[k]) for k in keys(bcontrols)];
        [(q₀[k], qcontrols[k]) for k in keys(qcontrols)];
        isnothing(mcontrols) ? [] : [(m₀, mcontrols)]]
    nb, nq, nm = map(length ∘ vec, values(u))
    ranges = (b = 1:nb, q = nb .+ (1:nq), m = nb + nq .+ (1:nm))

    ∇² = any(settings -> haskey(settings, :L), bcontrols) ? surfacelaplacianmatrix(γ) : nothing
    σ, lb, ub, smooth = Float64[], Float64[], Float64[], SparseMatrixCSC{Float64,Int}[]
    for (x₀, settings) in blocks
        x = vec(x₀)
        # a setting at each adjusted grid point: one number everywhere, or a field like x₀
        bygridpoint(v) = v isa Number ? fill(float(v), length(x)) : vec(v)
        σₖ = bygridpoint(settings.σ)
        logscale = get(settings, :logscale, false)
        # bound on the adjustment u that keeps x₀ + u (or x₀ eᵘ) within the
        # physical limit v: v - x₀, or log v - log x₀
        adjustmentbound(v) = logscale ? log.(v) - log.(x) : v - x
        append!(σ, σₖ)
        append!(lb, adjustmentbound(bygridpoint(get(settings, :lb, logscale ? 0.0 : -Inf))))
        append!(ub, adjustmentbound(bygridpoint(get(settings, :ub, Inf))))
        # smoothness of a boundary adjustment with length scale L: ∇²ᵀ Diagonal(L⁴/σ²) ∇²
        Q⁻ₖ = haskey(settings, :L) ?
            transpose(∇²) * spdiagm(0 => bygridpoint(settings.L) .^ 4 ./ σₖ .^ 2) * ∇² :
            spzeros(length(x), length(x))
        push!(smooth, (Q⁻ₖ + transpose(Q⁻ₖ)) / 2)
    end
    Q⁻ = spdiagm(0 => σ .^ -2) + blockdiag(smooth...)
    return Inversion(y, b₀, q₀, r, m₀, u, ranges, σ, lb, ub, Q⁻)
end

function Base.show(io::IO, inversion::Inversion)
    print(io, "Inversion($(length(inversion.σ)) controls, ",
        "$(length(inversion.y.tracers)) observed tracers)")
end

"""
    steadyinversion(Alu, b, q, r, γ)

Steady-state tracers for every boundary condition in `b`, one tracer per
boundary condition (`main`'s `steadyinversion(Alu, b::NamedTuple, γ; q, r)`
instead combines several boundaries of one tracer). The interior source of each
tracer is the sum of the sources in `q` scaled by its stoichiometric ratios in
`r`.

# Arguments
- `Alu`: LU factorization of the water-mass matrix
- `b::NamedTuple`: surface boundary conditions
- `q::NamedTuple`: interior sources
- `r::NamedTuple`: stoichiometric ratios, e.g. `(NO₃ = (qPO₄ = -15.5,),)`; a
  tracer missing from `r` has no interior source
- `γ::Grid`: TMI grid

# Output
- `c::NamedTuple`: steady-state `Field`s keyed like `b`
"""
function steadyinversion(Alu, b::NamedTuple, q::NamedTuple, r::NamedTuple, γ::Grid)
    return foreachtracer(keys(b)) do k
        rₖ = get(r, k, nothing)
        qₖ = isnothing(rₖ) ? nothing : reduce(+, map(j -> rₖ[j] * q[j], keys(rₖ)))
        steadyinversion(Alu, b[k], γ; q=qₖ)
    end
end

"""
    adjustfirstguess!(u, inversion, uvec)

The first-guess boundary conditions, sources, and mass fractions adjusted by
the control vector, which is written into the preallocated controls `u`.

# Arguments
- `u`: preallocated controls, a copy of `inversion.u`
- `inversion::Inversion`: first guesses and blocks of `uvec`
- `uvec`: control vector

# Output
- `(b, q, m)`: adjusted boundary conditions, sources, and mass fractions
"""
function adjustfirstguess!(u, inversion::Inversion, uvec::AbstractVector)
    ranges = inversion.ranges
    unvec!(u.b, uvec[ranges.b])
    unvec!(u.q, uvec[ranges.q])
    b = adjustboundarycondition(inversion.b₀, u.b)
    q = foreachtracer(keys(inversion.q₀)) do k
        haskey(u.q, k) ? adjustsource(inversion.q₀[k], u.q[k]) : inversion.q₀[k]
    end
    isempty(ranges.m) && return b, q, inversion.m₀
    # u.m holds the adjusted mass fractions, not an adjustment
    unvec!(u.m, vec(inversion.m₀) + uvec[ranges.m])
    return b, q, u.m
end

"""
    steadyinversion(uvec, inversion, u, F)

Steady-state tracers for the adjusted boundary conditions, sources, and mass
fractions of a control vector.

# Arguments
- `uvec`: control vector
- `inversion::Inversion`: first guesses
- `u`: preallocated controls, a copy of `inversion.u`
- `F`: UMFPACK factorization of a water-mass matrix on the same grid, refactored
  in place by `lu!`, which keeps its ordering

# Output
- `c::NamedTuple`: steady-state `Field`s keyed like `inversion.b₀`
"""
function steadyinversion(uvec::AbstractVector, inversion::Inversion, u, F)
    γ = inversion.y.γ
    b, q, m = adjustfirstguess!(u, inversion, uvec)
    return steadyinversion(lu!(F, watermassmatrix(m, γ)), b, q, inversion.r, γ)
end

"""
    costfunction(uvec, inversion, u, F)
    costterms(c, uvec, inversion)

Cost function `J = Jdata + Jcontrol`, the weighted squared model-data misfit
`Jdata = nᵀ Wⁱ n` for `n = ỹ - y` plus the control penalty
`Jcontrol = uvecᵀ Q⁻ uvec`.

# Arguments
- `uvec`: control vector
- `inversion::Inversion`: observations, first guesses, and weighting matrices
- `u`: preallocated controls, a copy of `inversion.u`; Enzyme differentiates
  through it
- `F`: UMFPACK factorization refactored in place, see `steadyinversion`
- `c::NamedTuple`: steady-state tracers for `uvec`

# Output
- `J`: cost function (`costfunction`)
- `(Jdata, Jcontrol)`: its two terms (`costterms`)
"""
function costterms(c::NamedTuple, uvec::AbstractVector, inversion::Inversion)
    y = inversion.y
    n = observe(c, y) - vec(y)
    return (Jdata = dot(n, y.Wⁱ, n), Jcontrol = dot(uvec, inversion.Q⁻, uvec))
end
costfunction(uvec::AbstractVector, inversion::Inversion, u, F) =
    sum(costterms(steadyinversion(uvec, inversion, u, F), uvec, inversion))

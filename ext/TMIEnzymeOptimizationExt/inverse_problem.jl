# The `Inversion` type: its constructor, how it prints, and the first guesses
# adjusted by a control vector. The struct itself is defined in src/inversion.jl.

"""
    function Inversion(y; b₀, q₀=(;), r=(;), m₀, controls)

Inverse problem for adjustments `u` to first-guess boundary conditions `b₀`,
sources `q₀`, and mass fractions `m₀` that fit observations `y`, with cost
`J = Jdata + Jb + Js + Jq + Jm` (see `costterms`).

# Arguments
- `y::Observations`: observations and their weighting matrix `Wⁱ`
- `b₀::NamedTuple`: first-guess surface `BoundaryCondition` for each tracer
- `m₀::NamedTuple`: first-guess `MassFraction`s
- `controls::NamedTuple`: adjusted quantities only, e.g.
  `(b = (θ = (σ=…, L=…, lb=…, ub=…),), q = (qPO₄ = (σ=…, lb=…, ub=…, logscale=true),), m = (σ=…, lb=…, ub=…))`.
  `σ` is the standard deviation of `u`, `L` [km] the smoothing scale of a
  boundary adjustment, and `lb`, `ub` are physical limits on `b`, `q`, and `m`;
  each is one number for every adjusted grid point or a field shaped like the
  first guess. A source with `logscale=true` is adjusted multiplicatively,
  `q = q₀ eᵘ`, so its `σ` is in log units (`σ = log(50)` is a one-σ factor of
  50), while its `lb` and `ub` remain limits on `q` itself.
# Optional Arguments
- `q₀::NamedTuple`: first-guess interior `Source`s, linear scale
- `r::NamedTuple`: stoichiometric ratios of each tracer to the sources, signed as
  `r` in `steadyinversion`, e.g. `(PO₄ = (qPO₄ = 1.0,), NO₃ = (qPO₄ = 15.5,))`
# Output
- `Inversion`: first guesses, preallocated controls `u`, blocks `ranges` of
  `uvec`, and `σ`, `lb`, `ub`, and the weights `smoothness` of the prior smoothness cost
"""
function Inversion(y::Observations; b₀::NamedTuple, q₀::NamedTuple=(;),
    r::NamedTuple=(;), m₀::NamedTuple, controls::NamedTuple)
    γ = y.γ
    bcontrols = get(controls, :b, (;))
    qcontrols = get(controls, :q, (;))
    mcontrols = get(controls, :m, nothing)

    # zero adjustments u for the controlled boundary conditions and sources; u.m
    # holds the adjusted mass fractions themselves, which watermassmatrix takes
    # (TMI has no arithmetic on MassFraction, so they are adjusted as a vector)
    ub₀ = map(b -> zerosurfaceboundary(γ, b.name, b.longname, b.units), b₀[keys(bcontrols)])
    uq₀ = map(q₀[keys(qcontrols)], qcontrols) do q, settings
        zerosource(γ, q.name, q.longname, q.units; logscale=get(settings, :logscale, false))
    end
    u = (b = ub₀, q = uq₀, m = isnothing(mcontrols) ? (;) : unvec(m₀, vec(m₀)))

    # controlled blocks of uvec, in b, q, m order
    blocks = [[(b₀[k], bcontrols[k]) for k in keys(bcontrols)];
        [(q₀[k], qcontrols[k]) for k in keys(qcontrols)];
        isnothing(mcontrols) ? [] : [(m₀, mcontrols)]]
    nb, nq, nm = map(uₖ -> length(vec(uₖ)), values(u))
    ranges = (b = 1:nb, q = nb .+ (1:nq), m = nb + nq .+ (1:nm))

    ∇² = any(settings -> haskey(settings, :L), bcontrols) ? surfacelaplacianmatrix(γ) : nothing
    σ, lb, ub, smooth = Float64[], Float64[], Float64[], SparseMatrixCSC{Float64,Int}[]
    for (x₀, settings) in blocks
        x = vec(x₀)
        logscale = get(settings, :logscale, false)
        # each setting at every adjusted point: one number for all of them, or a
        # field shaped like the first guess. Physical limits default to 0 for a
        # log-scale source and to unbounded otherwise.
        pointwise(v) = v isa Number ? fill(float(v), length(x)) : vec(v)
        σₖ = pointwise(settings.σ)
        lower = pointwise(get(settings, :lb, logscale ? 0.0 : -Inf))
        upper = pointwise(get(settings, :ub, Inf))
        Lₖ = haskey(settings, :L) ? pointwise(settings.L) : nothing

        # limits on the adjustment u that keep x₀ + u, or x₀ eᵘ, within the physical limits
        append!(σ, σₖ)
        if logscale
            append!(lb, log.(lower) - log.(x))
            append!(ub, log.(upper) - log.(x))
        else
            append!(lb, lower - x)
            append!(ub, upper - x)
        end

        # prior smoothness cost of a boundary adjustment with length scale L: ∇²ᵀ Diagonal(L⁴/σ²) ∇²
        if isnothing(Lₖ)
            smoothₖ = spzeros(length(x), length(x))
        else
            smoothₖ = transpose(∇²) * spdiagm(0 => Lₖ .^ 4 ./ σₖ .^ 2) * ∇²
        end
        # symmetric in exact arithmetic; averaging with its transpose removes round-off
        push!(smooth, (smoothₖ + transpose(smoothₖ)) / 2)
    end
    return Inversion(y, b₀, q₀, r, m₀, u, ranges, σ, lb, ub, blockdiag(smooth...))
end

function Base.show(io::IO, inversion::Inversion)
    print(io, "Inversion($(length(inversion.σ)) controls, ",
        "$(length(inversion.y.values)) observed tracers)")
end

"""
    function adjustfirstguess!(u, inversion, uvec)

The first-guess boundary conditions, sources, and mass fractions adjusted by
the control vector, which is written into the preallocated controls `u`. Kept
apart from `forwardmodel`: written inside it, Enzyme's type analysis fails
(`IllegalTypeAnalysisException`) for inversions with fixed tracers or sources.

# Arguments
- `u`: preallocated controls, a copy of `inversion.u`: the adjustments `u.b` and
  `u.q` as `BoundaryCondition`s and `Source`s, and the adjusted mass fractions `u.m`
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
    # TMI's adjustsource on a NamedTuple needs an adjustment for every source, so
    # only the adjusted sources go through it and the fixed ones are kept as given
    q = merge(inversion.q₀, map(adjustsource, inversion.q₀[keys(u.q)], u.q))
    isempty(ranges.m) && return b, q, inversion.m₀
    # u.m holds the adjusted mass fractions, not an adjustment: watermassmatrix
    # takes the fractions, and writing them into u.m allocates no new MassFractions
    unvec!(u.m, vec(inversion.m₀) + uvec[ranges.m])
    return b, q, u.m
end

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
    u = (b = foreachtracer(k -> zerosurfaceboundary(γ, b₀[k].name, b₀[k].longname, b₀[k].units),
            keys(bcontrols)),
        q = foreachtracer(k -> zerosource(γ, q₀[k].name, q₀[k].longname, q₀[k].units;
            logscale=get(qcontrols[k], :logscale, false)), keys(qcontrols)),
        m = isnothing(mcontrols) ? (;) : unvec(m₀, vec(m₀)))

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
        # a setting at each adjusted grid point: one number everywhere, or a field like x₀
        bygridpoint(v) = v isa Number ? fill(float(v), length(x)) : vec(v)
        σₖ = bygridpoint(settings.σ)
        logscale = get(settings, :logscale, false)
        # bound on the adjustment u that keeps x₀ + u (or x₀ eᵘ) within the
        # physical limit v: v - x₀, or log v - log x₀
        adjustmentbound(v) = logscale ? log.(v) - log.(x) : v - x
        # physical limits, by default 0 for a log-scale source and unbounded otherwise
        lower = get(settings, :lb, logscale ? 0.0 : -Inf)
        upper = get(settings, :ub, Inf)
        append!(σ, σₖ)
        append!(lb, adjustmentbound(bygridpoint(lower)))
        append!(ub, adjustmentbound(bygridpoint(upper)))
        # prior smoothness cost of a boundary adjustment with length scale L: ∇²ᵀ Diagonal(L⁴/σ²) ∇²
        smoothₖ = haskey(settings, :L) ?
            transpose(∇²) * spdiagm(0 => bygridpoint(settings.L) .^ 4 ./ σₖ .^ 2) * ∇² :
            spzeros(length(x), length(x))
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
    q = foreachtracer(keys(inversion.q₀)) do k
        haskey(u.q, k) ? adjustsource(inversion.q₀[k], u.q[k]) : inversion.q₀[k]
    end
    isempty(ranges.m) && return b, q, inversion.m₀
    # u.m holds the adjusted mass fractions, not an adjustment: watermassmatrix
    # takes the fractions, and writing them into u.m allocates no new MassFractions
    unvec!(u.m, vec(inversion.m₀) + uvec[ranges.m])
    return b, q, u.m
end

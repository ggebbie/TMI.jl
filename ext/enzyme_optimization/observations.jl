# an observation setting (σ, locs, L, Wⁱ) for each tracer in `names`: one value
# for all tracers, or a `NamedTuple` by tracer (`nothing` for tracers missing from it)
bytracer(a, names) = foreachtracer(k -> a isa NamedTuple ? get(a, k, nothing) : a, names)

"""
    atobservations(c, wis, γ, n)

Values of a tracer, or of its σ or L, at the points where the tracer is observed.

# Arguments
- `c`: gridded `Field`, vector of observed values, or number
- `wis`: interpolation weights at observation locations, or `nothing` for gridded observations
- `γ::Grid`: TMI grid
- `n`: number of observations, used when `c` is a number

# Output
- `ỹ`: interior cells of a gridded `Field`, `observe(c, wis, γ)` at observation
  locations, a vector as given, or a number repeated `n` times
"""
atobservations(c::Field, ::Nothing, γ, n=0) = c.tracer[γ.interior]
atobservations(c::Field, wis, γ, n=0) = observe(c, wis, γ)
atobservations(c::AbstractVector, wis, γ, n=0) = c
atobservations(c::Number, wis, γ, n) = fill(float(c), n)

"""
    observe(c, y)

Values of the tracers `c` where `y` observes them.

# Arguments
- `c::NamedTuple`: gridded `Field`s, including every tracer observed in `y`
- `y::Observations`: observations

# Output
- `ỹ`: vector ordered like `vec(y)`
"""
observe(c::NamedTuple, y::Observations) =
    reduce(vcat, map(k -> atobservations(c[k], y.wis[k], y.γ), keys(y.tracers)))

"""
    Observations(y, σ, γ; locs=nothing, L=nothing, Wⁱ=nothing)

Tracer observations and their error statistics. `σ`, `locs`, `L`, and `Wⁱ` are
given once for all tracers or as a `NamedTuple` keyed like `y`.

# Arguments
- `y::NamedTuple`: each tracer a gridded `Field` or a vector of values at `locs`
- `σ`: observational error, a number, `Field`, or vector
- `γ::Grid`: TMI grid
- `locs`: `nothing` for gridded tracers, or a vector of (lon, lat, depth)
- `L`: decorrelation length scale [km] of gridded errors, a number or `Field`
- `Wⁱ`: inverse observational error covariance, overrides `σ` and `L`

# Output
- `Observations`: `tracers` (the given `y`), `σ`, `locs`, interpolation weights
  `wis`, all observed values `yvec` (also `vec(y)`), `γ`, and block-diagonal `Wⁱ`
"""
function Observations(y::NamedTuple, σ, γ::Grid; locs=nothing, L=nothing, Wⁱ=nothing)
    names = keys(y)
    σ, locs, L, Wⁱ = map(a -> bytracer(a, names), (σ, locs, L, Wⁱ))
    wis = map(x -> isnothing(x) ? nothing : [TMI.interpindex(loc, γ) for loc in x], locs)
    observed = map(k -> atobservations(y[k], wis[k], γ), names)
    Wⁱ = map(names, observed) do k, yₖ
        n = length(yₖ)
        if !isnothing(Wⁱ[k])
            sparse(Wⁱ[k])
        elseif !isnothing(L[k])
            isnothing(wis[k]) || throw(ArgumentError(
                "L needs gridded observations of $k; supply Wⁱ for observations at locs"))
            gaussianprecision(atobservations(σ[k], nothing, γ, n), atobservations(L[k], nothing, γ, n), γ)
        else
            spdiagm(0 => atobservations(σ[k], wis[k], γ, n) .^ -2)
        end
    end
    return Observations(y, σ, locs, wis, reduce(vcat, observed), γ, blockdiag(Wⁱ...))
end

Base.vec(y::Observations) = y.yvec

function Base.show(io::IO, y::Observations)
    gridded = count(isnothing, y.wis)
    print(io, "Observations($(length(y.tracers)) tracers; $gridded gridded, ",
        "$(length(y.tracers) - gridded) at locations; $(length(vec(y))) values)")
end

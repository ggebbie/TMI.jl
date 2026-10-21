# an observation setting (σ, locs, L, Wⁱ) for each tracer in `names`: one value
# for all tracers, or a `NamedTuple` by tracer (`nothing` for tracers missing from it)
bytracer(a, names) = foreachtracer(k -> a isa NamedTuple ? get(a, k, nothing) : a, names)

"""
    function observe(c, y)

Values of the tracers `c` where `y` observes them: at the interior cells for a
gridded tracer, interpolated to the observation locations otherwise.

# Arguments
- `c::NamedTuple`: gridded `Field`s, including every tracer observed in `y`
- `y::Observations`: observations
# Output
- `ỹ`: vector ordered like `vec(y)`
"""
observe(c::NamedTuple, y::Observations) = reduce(vcat, map(keys(y.values)) do k
    isnothing(y.wis[k]) ? c[k].tracer[y.γ.interior] : observe(c[k], y.wis[k], y.γ)
end)

"""
    function Observations(observed, σ, γ; locs=nothing, L=nothing, Wⁱ=nothing)

Tracer observations and their error statistics. `σ`, `locs`, `L`, and `Wⁱ` are
given once for all tracers or as a `NamedTuple` keyed like `observed`.

# Arguments
- `observed::NamedTuple`: the observed tracers, each a gridded `Field` or a
  vector of values at `locs`
- `σ`: observational error, a number, `Field`, or vector
- `γ::Grid`: TMI grid
# Optional Arguments
- `locs`: `nothing` for gridded tracers, or a vector of (lon, lat, depth)
- `L`: decorrelation length scale [km] of gridded errors, a number or `Field`
- `Wⁱ`: inverse observational error covariance, overrides `σ` and `L`
# Output
- `Observations`: for each tracer, its observed `values` and their `σ` at the
  observation points (the interior cells of a gridded tracer); `locs`,
  interpolation weights `wis`, all observed values `yvec` (also `vec(y)`), `γ`,
  and block-diagonal `Wⁱ`
"""
function Observations(observed::NamedTuple, σ, γ::Grid; locs=nothing, L=nothing, Wⁱ=nothing)
    names = keys(observed)
    σ, locs, L, Wⁱ = map(a -> bytracer(a, names), (σ, locs, L, Wⁱ))
    wis = map(x -> isnothing(x) ? nothing : [TMI.interpindex(loc, γ) for loc in x], locs)
    # a setting of tracer k at each of its n observation points: a Field at the
    # interior cells or interpolated to `locs`, a vector as given, a number repeated
    function atpoints(x, k, n=0)
        x isa Number && return fill(float(x), n)
        x isa Field || return x
        return isnothing(wis[k]) ? x.tracer[γ.interior] : observe(x, wis[k], γ)
    end
    pointvalues = foreachtracer(k -> atpoints(observed[k], k), names)
    pointσ = foreachtracer(k -> atpoints(σ[k], k, length(pointvalues[k])), names)
    Wⁱ = map(names) do k
        if !isnothing(Wⁱ[k])
            sparse(Wⁱ[k])
        elseif !isnothing(L[k])
            isnothing(wis[k]) || throw(ArgumentError("L applies to gridded observations, " *
                "and $k is observed at locs: give L by tracer, e.g. L = (θ = 1000.0,), or supply Wⁱ"))
            gaussianprecision(pointσ[k], atpoints(L[k], k, length(pointvalues[k])), γ)
        else
            spdiagm(0 => pointσ[k] .^ -2)
        end
    end
    return Observations(pointvalues, pointσ, locs, wis, reduce(vcat, pointvalues), γ, blockdiag(Wⁱ...))
end

Base.vec(y::Observations) = y.yvec

function Base.show(io::IO, y::Observations)
    gridded = count(isnothing, y.wis)
    print(io, "Observations($(length(y.values)) tracers; $gridded gridded, ",
        "$(length(y.values) - gridded) at locations; $(length(vec(y))) values)")
end

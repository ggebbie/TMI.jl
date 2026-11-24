# The `Observations` type: its constructor, sampling modeled tracers where they are
# observed, and how it prints. The struct itself is defined in src/inversion.jl.

"""
    function observe(c, y)

Method of TMI's `observe` for all tracers of an `Observations`. Values of the
tracers `c` where `y` observes them: at the interior cells for a gridded tracer,
interpolated to the observation locations otherwise.

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
    # each setting by tracer: one value for all tracers, or a NamedTuple by tracer
    # (`nothing` for tracers missing from it)
    σbytracer, locsbytracer, Lbytracer, Wⁱbytracer = map((σ, locs, L, Wⁱ)) do a
        foreachtracer(k -> a isa NamedTuple ? get(a, k, nothing) : a, names)
    end
    # each tracer's observation points: its locations, with their interpolation
    # weights, or the interior cells if it is gridded
    wis = map(x -> isnothing(x) ? nothing : [TMI.interpindex(loc, γ) for loc in x], locsbytracer)
    npoints = map(w -> isnothing(w) ? sum(γ.interior) : length(w), wis)
    # a setting of one tracer at each of its n observation points: a Field at the
    # interior cells or interpolated to its locations, a vector as given, a number repeated
    function atpoints(x, wisₖ, n)
        if x isa Number
            return fill(float(x), n)
        elseif x isa Field
            return isnothing(wisₖ) ? x.tracer[γ.interior] : observe(x, wisₖ, γ)
        else
            return x
        end
    end
    pointvalues = map(atpoints, observed, wis, npoints)
    pointσ = map(atpoints, σbytracer, wis, npoints)
    Wⁱblocks = map(names) do k
        if !isnothing(Wⁱbytracer[k])
            sparse(Wⁱbytracer[k])
        elseif !isnothing(Lbytracer[k])
            if !isnothing(wis[k])
                throw(ArgumentError("L applies to gridded observations, and $k is observed " *
                    "at locs: give L by tracer, e.g. L = (θ = 1000.0,), or supply Wⁱ"))
            end
            gaussianprecision(pointσ[k], atpoints(Lbytracer[k], wis[k], npoints[k]), γ)
        else
            spdiagm(0 => pointσ[k] .^ -2)
        end
    end
    return Observations(pointvalues, pointσ, locsbytracer, wis, reduce(vcat, pointvalues), γ,
        blockdiag(Wⁱblocks...))
end

Base.vec(y::Observations) = y.yvec

function Base.show(io::IO, y::Observations)
    gridded = count(isnothing, y.wis)
    print(io, "Observations($(length(y.values)) tracers; $gridded gridded, ",
        "$(length(y.values) - gridded) at locations; $(length(vec(y))) values)")
end

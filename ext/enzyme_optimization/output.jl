"""
    writeatomically(f, file)

Write a file through a temporary copy, so that an interrupted run never leaves a
partial file behind.

# Arguments
- `f`: function that writes to the file name it is given
- `file`: file name

# Output
- `file`
"""
function writeatomically(f, file)
    mkpath(dirname(file))
    tmp = file * ".tmp"
    rm(tmp; force=true)
    f(tmp)
    mv(tmp, file; force=true)
    return file
end

"""
    griddedobservations(values, locs, γ, name, longname, units)

Observations placed on the TMI grid: each value at its nearest grid cell.

# Arguments
- `values`: observed values
- `locs`: `nothing` for values at interior cells, or (lon, lat, depth) of each value
- `γ::Grid`: TMI grid
- `name`, `longname`, `units`: `Field` metadata

# Output
- `Field`: the values at interior cells, or the mean of the values nearest to
  each grid cell, with NaN elsewhere
"""
function griddedobservations(values, ::Nothing, γ, name, longname, units)
    tracer = fill(NaN, size(γ.wet))
    tracer[γ.interior] .= values
    return Field(tracer, γ, name, longname, units)
end

function griddedobservations(values, locs, γ, name, longname, units)
    total, count = zeros(size(γ.wet)), zeros(Int, size(γ.wet))
    for (value, loc) in zip(values, locs)
        cell = TMI.nearestneighbor(loc, γ)
        total[cell] += value
        count[cell] += 1
    end
    tracer = fill(NaN, size(γ.wet))
    observed = count .> 0
    tracer[observed] .= total[observed] ./ count[observed]
    return Field(tracer, γ, name, longname, units)
end

"""
    saveobservations(name, directory, inversion)

Write the observations to `directory/name_observations.nc`. Each tracer `k` and
its `σk` are written as gridded fields; `k` has the attribute `support`
(`"gridded"` or `"sparse"`). Tracers observed at locations also get their
values, `σ`, and coordinates as `k_sample`, `σk_sample`, `k_longitude`,
`k_latitude`, and `k_depth`.

# Arguments
- `name`: experiment name
- `directory`: output directory
- `inversion::Inversion`: holds the observations

# Output
- `file`: NetCDF file name
"""
function saveobservations(name, directory, inversion::Inversion)
    y = inversion.y
    γ = y.γ
    file = joinpath(directory, "$(name)_observations.nc")
    writeatomically(file) do tmp
        for k in keys(y.tracers)
            locs, wis = y.locs[k], y.wis[k]
            yₖ = atobservations(y.tracers[k], wis, γ)
            σₖ = atobservations(y.σ[k], wis, γ, length(yₖ))
            longname, units = inversion.b₀[k].longname, inversion.b₀[k].units
            σname = "σ$k"
            writefield(tmp, griddedobservations(yₖ, locs, γ, k, longname, units))
            writefield(tmp, griddedobservations(σₖ, locs, γ, Symbol(σname),
                "1σ uncertainty in $longname", units))
            NCDataset(tmp, "a") do ds
                ds[String(k)].attrib["support"] = isnothing(locs) ? "gridded" : "sparse"
                isnothing(locs) && return
                samples = "$(k)_sample"
                defDim(ds, samples, length(yₖ))
                for (vname, values, vlongname, vunits) in ((samples, yₖ, longname, units),
                    ("$(σname)_sample", σₖ, "1σ uncertainty in $longname", units),
                    ("$(k)_longitude", first.(locs), "longitude", "°E"),
                    ("$(k)_latitude", getindex.(locs, 2), "latitude", "°N"),
                    ("$(k)_depth", last.(locs), "depth", "m"))
                    defVar(ds, vname, Float64.(values), (samples,);
                        attrib=Dict("longname" => vlongname, "units" => vunits))
                end
            end
        end
    end
    println("Saved observations to $file")
    return file
end

"""
    saveinversion(name, directory, uvec, cache; iteration, J, final=false)

Write the modeled tracers, sources, water-mass matrix, and cost function for a
control vector to NetCDF.

# Arguments
- `name`: experiment name
- `directory`: output directory
- `uvec`: control vector
- `cache::EnzymeCostGradientCache`: the inverse problem, and the controls `u`
  and factorization `F` that are overwritten here
- `iteration`: optimizer iteration
- `J`: cost function at `uvec`
- `final`: name the file `name_final.nc` instead of `name_iteration_NNNNNN.nc`

# Output
- `file`: NetCDF file name
"""
function saveinversion(name, directory, uvec, cache;
    iteration::Integer, J::Real, final::Bool=false)
    (; inversion, u, F) = cache
    γ = inversion.y.γ
    b, q, m = adjustfirstguess!(u, inversion, uvec)
    A = watermassmatrix(m, γ)
    c = steadyinversion(lu!(F, A), b, q, inversion.r, γ)
    terms = costterms(c, uvec, inversion)
    suffix = final ? "final" : "iteration_$(lpad(iteration, 6, '0'))"
    file = joinpath(directory, "$(name)_$suffix.nc")
    writeatomically(file) do tmp
        foreach(k -> writefield(tmp, Field(c[k].tracer, γ, k, c[k].longname, c[k].units)), keys(c))
        foreach(qₖ -> TMI.writesource(tmp, qₖ), q)
        TMI.watermassmatrix2nc("", A; filenetcdf=tmp)
        NCDataset(tmp, "a") do ds
            ds.attrib["experiment_name"] = name
            ds.attrib["iteration"] = iteration
            ds.attrib["J"] = J
            ds.attrib["Jdata"] = terms.Jdata
            ds.attrib["Jcontrol"] = terms.Jcontrol
        end
    end
    println("Saved iteration $iteration to $file")
    return file
end

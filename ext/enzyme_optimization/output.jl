"""
    function save(file, inversion)
    function save(file, workspace, uvec; iteration)

Write to the NetCDF `file` either the observations an inversion fits, or the
state of the inversion at the control vector `uvec`.
- `save(file, inversion)`: for each tracer `k`, its observed values `k` and
  their uncertainty `σk`, at the interior cells (in `γ.interior` order) for
  gridded tracers, and with the coordinates `k_longitude`, `k_latitude`, and
  `k_depth` for tracers observed at locations.
- `save(file, workspace, uvec; iteration)`: the modeled tracers, the sources, the
  water-mass matrix, and the attributes `iteration`, `J`, its data cost `Jdata` and
  prior costs `Jb`, `Js`, `Jq`, and `Jm` (see `costterms`), and their sum
  `Jcontrol = Jb + Js + Jq + Jm`.

# Arguments
- `file`: NetCDF file name, replaced if it exists
- `inversion::Inversion`: holds the observations
- `workspace::InversionWorkspace`: the inverse problem, and the controls `u` and
  factorization `F` that are overwritten here
- `uvec`: control vector
- `iteration`: optimizer iteration
# Output
- `file`
"""
function save(file, inversion::Inversion)
    y = inversion.y
    mkpath(dirname(file))
    rm(file; force=true)
    NCDataset(file, "c") do ds
        for k in keys(y.values)
            locs = y.locs[k]
            longname, units = inversion.b₀[k].longname, inversion.b₀[k].units
            variables = [(String(k), y.values[k], longname, units),
                ("σ$k", y.σ[k], "1σ uncertainty in $longname", units)]
            isnothing(locs) || append!(variables, [
                ("$(k)_longitude", first.(locs), "longitude", "°E"),
                ("$(k)_latitude", getindex.(locs, 2), "latitude", "°N"),
                ("$(k)_depth", last.(locs), "depth", "m")])
            observation = "$(k)_observation"
            defDim(ds, observation, length(y.values[k]))
            for (vname, data, vlongname, vunits) in variables
                defVar(ds, vname, Float64.(data), (observation,);
                    attrib=Dict("longname" => vlongname, "units" => vunits))
            end
        end
    end
    println("Saved observations to $file")
    return file
end

function save(file, workspace::InversionWorkspace, uvec; iteration::Integer)
    (; inversion, u, F) = workspace
    γ = inversion.y.γ
    (; q, A, c) = forwardmodel(uvec, inversion, u, F)
    terms = costterms(c, uvec, inversion)
    mkpath(dirname(file))
    rm(file; force=true)
    foreach(k -> writefield(file, Field(c[k].tracer, γ, k, c[k].longname, c[k].units)), keys(c))
    foreach(qₖ -> TMI.writesource(file, qₖ), q)
    TMI.watermassmatrix2nc("", A; filenetcdf=file)
    NCDataset(file, "a") do ds
        ds.attrib["iteration"] = iteration
        ds.attrib["J"] = sum(terms)
        foreach(name -> ds.attrib[String(name)] = terms[name], keys(terms))
        ds.attrib["Jcontrol"] = terms.Jb + terms.Js + terms.Jq + terms.Jm
    end
    println("Saved iteration $iteration to $file")
    return file
end

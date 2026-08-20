"""
    InversionCheckpointer(name, interval; checkpoint_directory, data_directory)

Configure inversion checkpoint filenames and frequency.

# Arguments
- `name`, `interval`, `checkpoint_directory`, `data_directory`: output configuration

# Output
- `checkpointer`: checkpoint configuration
"""
struct InversionCheckpointer
    name::String
    interval::Int
    checkpoint_directory::String
    data_directory::String
end
function InversionCheckpointer(name::AbstractString, interval::Integer;
    checkpoint_directory::AbstractString=TMI.pkgdir("checkpoints", name),
    data_directory::AbstractString=TMI.pkgdatadir())
    return InversionCheckpointer(String(name), interval,
        String(checkpoint_directory), String(data_directory))
end

"""
    sparse2field(y, locations, γ, name, longname, units)

Convert sparse or gridded observation data to a diagnostic field.

# Arguments
- `y`, `locations`, `γ`: data, optional coordinates, and target grid
- `name`, `longname`, `units`: output field metadata

# Output
- `field`: diagnostic field, or `nothing` when `y` is absent
"""
function sparse2field(y::AbstractVector, locations::AbstractVector,
    γ::Grid, name, longname, units)
    tracer = fill(NaN, size(γ.wet))
    sample_sum = zeros(size(tracer))
    sample_count = zeros(Int, size(tracer))
    for (value, location) in zip(y, locations)
        cell = TMI.nearestneighbor(location, γ)
        sample_sum[cell] += value
        sample_count[cell] += 1
    end
    sampled_cells = sample_count .> 0
    tracer[sampled_cells] .=
        sample_sum[sampled_cells] ./ sample_count[sampled_cells]
    return Field(tracer, γ, name, longname, units)
end
function sparse2field(y::Field, locations::AbstractVector,
    γ::Grid, name, longname, units)
    return sparse2field(observe(y, interpindex(locations, γ), γ), locations,
        γ, name, longname, units)
end
function sparse2field(y::Real, locations::AbstractVector,
    γ::Grid, name, longname, units)
    return sparse2field(fill(y, length(locations)), locations,
        γ, name, longname, units)
end
function sparse2field(y::Field, ::Nothing, ::Grid, name, longname, units)
    return Field(y.tracer, y.γ, name, longname, units)
end
function sparse2field(y::AbstractVector, ::Nothing,
    γ::Grid, name, longname, units)
    tracer = fill(NaN, size(γ.wet))
    tracer[γ.interior] .= y
    return Field(tracer, γ, name, longname, units)
end
sparse2field(y::Real, ::Nothing, γ::Grid, name, longname, units) =
    sparse2field(fill(y, sum(γ.interior)), nothing, γ, name, longname, units)

"""
    saveobservations(checkpointer, inversion)

Write observation fields, sparse samples, coordinates, and metadata to NetCDF.

# Arguments
- `checkpointer`, `inversion`: output configuration and inversion description

# Output
- `filename`: observation NetCDF filename
"""
function saveobservations(checkpointer::InversionCheckpointer, inversion::Inversion)
    observations = inversion.observations
    γ = observations.γ
    filename = joinpath(checkpointer.checkpoint_directory, "$(checkpointer.name)_observations.nc")
    temporary_filename = filename * ".tmp"
    mkpath(dirname(filename))
    rm(temporary_filename; force=true)
    writefield(temporary_filename,
        Field(γ.wet, γ, :wet, "wet grid cell", "1"))
    written_variable_names = Set(["wet"])
    NCDataset(dataset -> begin
        dataset.attrib["experiment_name"] = checkpointer.name
        dataset.attrib["title"] = "TMI inversion observations"
    end, temporary_filename, "a")

    for tracer_name in keys(observations.y)
        tracer_y = observations.y[tracer_name]
        tracer_σ = observations.σ[tracer_name]
        tracer_L = observations.L[tracer_name]
        locations = observations.locations[tracer_name]
        interpolation_indices = observations.interpolation_indices[tracer_name]
        tracer_template = inversion.control_templates.boundary_conditions[tracer_name]

        y_field = sparse2field(tracer_y, locations, γ, tracer_template.name,
            tracer_template.longname, tracer_template.units)
        σ_name = tracer_σ isa Field ? tracer_σ.name :
            Symbol("σ", tracer_template.name)
        σ_field = sparse2field(tracer_σ, locations, γ, σ_name,
            "1σ uncertainty in $(tracer_template.longname)", tracer_template.units)
        if isnothing(tracer_L)
            L_field = nothing
        else
            L_name = tracer_L isa Field ? tracer_L.name :
                Symbol("L", tracer_template.name)
            L_field = sparse2field(tracer_L, locations, γ, L_name,
                "decorrelation scale for $(tracer_template.longname)", "km")
        end

        diagnostic_variable_names = String[]
        for diagnostic_field in (y_field, σ_field, L_field)
            isnothing(diagnostic_field) && continue
            variable_name = String(diagnostic_field.name)
            if !(variable_name in written_variable_names)
                writefield(temporary_filename, diagnostic_field)
            end
            push!(written_variable_names, variable_name)
            push!(diagnostic_variable_names, variable_name)
        end
        y_variable_name, σ_variable_name = diagnostic_variable_names[1:2]
        L_variable_name = isnothing(L_field) ? nothing : String(L_field.name)

        NCDataset(dataset -> begin
            y_variable = dataset[y_variable_name]
            y_variable.attrib["model_variable"] = String(tracer_template.name)
            y_variable.attrib["uncertainty"] = σ_variable_name
            if isnothing(locations)
                y_variable.attrib["support"] = "gridded"
                if !isnothing(L_variable_name)
                    y_variable.attrib["decorrelation_scale"] = L_variable_name
                end
            else
                y_variable.attrib["support"] = "sparse"
                y_values = tracer_y isa Field ?
                    observe(tracer_y, interpolation_indices, γ) : tracer_y
                if tracer_σ isa Real
                    σ_values = fill(tracer_σ, length(y_values))
                elseif tracer_σ isa Field
                    σ_values = observe(tracer_σ, interpolation_indices, γ)
                else
                    σ_values = tracer_σ
                end
                if isnothing(tracer_L)
                    L_values = nothing
                elseif tracer_L isa Real
                    L_values = fill(tracer_L, length(y_values))
                elseif tracer_L isa Field
                    L_values = observe(tracer_L, interpolation_indices, γ)
                else
                    L_values = tracer_L
                end
                sample_dimension = "$(y_variable_name)_sample"
                defDim(dataset, sample_dimension, length(y_values))
                sample_variables = (
                    (sample_dimension, y_values, tracer_template.longname,
                        tracer_template.units),
                    ("$(σ_variable_name)_sample", σ_values,
                        "1σ uncertainty in $(tracer_template.longname)",
                        tracer_template.units),
                    ("$(y_variable_name)_longitude", first.(locations),
                        "longitude", "°E"),
                    ("$(y_variable_name)_latitude", getindex.(locations, 2),
                        "latitude", "°N"),
                    ("$(y_variable_name)_depth", last.(locations), "depth", "m"),
                )
                for (name, values, description, units) in sample_variables
                    sample_variable = defVar(dataset, name, Float64,
                        (sample_dimension,); attrib=Dict(
                            "longname" => description, "units" => units))
                    sample_variable[:] = values
                end
                y_variable.attrib["samples"] = sample_dimension
                if !isnothing(L_values)
                    scale_variable = defVar(dataset,
                        "$(L_variable_name)_sample", Float64,
                        (sample_dimension,); attrib=Dict(
                            "longname" => "decorrelation scale", "units" => "km"))
                    scale_variable[:] = L_values
                    y_variable.attrib["decorrelation_scale"] = L_variable_name
                end
            end
        end, temporary_filename, "a")
    end
    mv(temporary_filename, filename; force=true)
    println("Saved observations to $filename")
    return filename
end

"""
    saveinversion(checkpointer, x, inversion; iteration, objective, terminal=false)

Write modeled tracers, sources, transport, and optimizer state to NetCDF.

# Arguments
- `checkpointer`, `x`, `inversion`: output configuration, controls, and inversion
- `iteration`, `objective`, `terminal`: optimizer state and final-file selector

# Output
- `filename`: checkpoint or final NetCDF filename
"""
function saveinversion(checkpointer::InversionCheckpointer, x, inversion::Inversion;
    iteration::Integer, objective::Real, terminal::Bool=false)
    controls = unvec(inversion, x)
    γ = inversion.observations.γ
    A = watermassmatrix(controls.mass_fractions, γ)
    modeled_fields = steadyinversion(lu(A), controls.boundary_conditions,
        controls.sources, inversion.stoichiometry, γ)
    directory = terminal ? checkpointer.data_directory :
        checkpointer.checkpoint_directory
    filename_suffix = terminal ? "final" :
        "iteration_$(lpad(iteration, 6, '0'))"
    filename = joinpath(directory, "$(checkpointer.name)_$filename_suffix.nc")
    mkpath(dirname(filename))
    temporary_filename = filename * ".tmp"
    rm(temporary_filename; force=true)
    foreach(field -> writefield(temporary_filename, field), modeled_fields)
    foreach(source -> TMI.writesource(temporary_filename, source), controls.sources)
    TMI.watermassmatrix2nc("", A; file=temporary_filename)
    NCDataset(dataset -> begin
        dataset.attrib["experiment_name"] = checkpointer.name
        dataset.attrib["iteration"] = iteration
        dataset.attrib["objective"] = objective
    end, temporary_filename, "a")
    mv(temporary_filename, filename; force=true)
    println("Saved inversion iteration $iteration to $filename")
    return filename
end

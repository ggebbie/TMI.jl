#= Plot a saved inversion with `JULIA_LOAD_PATH="@v1.11:.:@stdlib" julia --project=@v1.11 scripts/inversion_diagnostics.jl full_inversion`. =#
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = get(ENV, "PYTHON", "python")
ENV["MPLBACKEND"] = "Agg"
using GeoPythonPlot, LinearAlgebra, NCDatasets, Statistics, TMI
const plt = GeoPythonPlot.pyplot
experiment = isempty(ARGS) ? "full_inversion" : only(ARGS)
reference = "modern_90x45x33_G14_v2"
checkpointdir, plotdir = TMI.pkgdir("checkpoints", experiment), TMI.pkgplotsdir(experiment)
files = (observations=joinpath(checkpointdir, "$(experiment)_observations.nc"), initial=joinpath(checkpointdir, "$(experiment)_iteration_000000.nc"),
    final=TMI.pkgdatadir("$(experiment)_final.nc"), reference=TMI.pkgdatadir("TMI_$reference.nc"))
mkpath(plotdir); println("Saving figures to $plotdir")
γ = TMI.Grid(files.observations, "wet", "lon", "lat", "depth"; flipdepth=false)
final_iteration = NCDataset(ds -> Int(ds.attrib["iteration"]), files.final)
"""
    plotfile(name); plottracer(tracer); volumemap(axis, field, title)

Construct output paths and the tracer and volume diagnostic panels.

# Arguments
- `name`, `tracer`, `axis`, `field`, `title`: output and plotting inputs

# Output
- `result`: file path, saved tracer figure, or plotted volume map
"""
plotfile(name) = joinpath(plotdir, "$(experiment)_$name.png")
# Observed, initial, and final tracer fields.
tracers = ((name="θ", label="Temperature", units="°C"), (name="Sₚ", label="Salinity", units="PSS-78"),
    (name="O₂", label="Oxygen", units="μmol/kg"))
function plottracer(tracer)
    fields = (readfield(files.observations, tracer.name, γ), readfield(files.initial, tracer.name, γ), readfield(files.final, tracer.name, γ))
    surface = argmin(abs.(γ.depth))
    rows = (("Surface", γ.lon, γ.lat, map(f -> f.tracer[:, :, surface], fields)), ("Zonal mean", γ.lat, -γ.depth, map(TMI.zonalaverage, fields)))
    titles = ("Observations", "Iteration 0", "Iteration $final_iteration", "Final − iteration 0")
    fig, ax = plt.subplots(2; ncols=4, figsize=(22, 9), sharex="row", sharey="row", squeeze=false, constrained_layout=true)
    for (row, (view, x, y, states)) in enumerate(rows)
        limits = extrema(filter(isfinite, vcat(states...)))
        difference = states[3] .- states[2]
        difference_limit = max(maximum(abs, filter(isfinite, difference)), eps())
        maps = []
        for (column, (title, values)) in enumerate(zip(titles, (states..., difference)))
            isdifference = column == 4
            levels = isdifference ? range(-difference_limit, difference_limit; length=21) : range(limits...; length=21)
            axᵢ = ax[row - 1, column - 1]
            push!(maps, axᵢ.contourf(x, y, values'; levels, cmap=isdifference ? "RdBu_r" : nothing))
            axᵢ.set(title="$(tracer.label) $view — $title", xlabel=view == "Surface" ? "Longitude" : "Latitude",
                ylabel=column == 1 ? (view == "Surface" ? "Latitude" : "Depth [m]") : "")
        end
        fig.colorbar(maps[3]; ax=[ax[row - 1, column] for column in 0:2], label=tracer.units)
        fig.colorbar(maps[4]; ax=ax[row - 1, 3], label="$(tracer.units) difference")
    end
    fig.savefig(plotfile("$(tracer.name)_quick_diagnostics"); dpi=200); plt.close(fig)
end
foreach(plottracer, tracers)
# Objective and normalized tracer misfits.
observation_specs = NCDataset(ds -> begin
    names = filter(name -> haskey(ds[name].attrib, "support"), keys(ds))
    map(function (name)
        variable, σname = ds[name], String(ds[name].attrib["uncertainty"])
        model = String(variable.attrib["model_variable"])
        if variable.attrib["support"] == "gridded"
            return (; name, model, values=variable[:, :, :][γ.interior], σ=ds[σname][:, :, :][γ.interior], locations=nothing)
        end
        sample = String(variable.attrib["samples"])
        locations = collect(zip(ds["$(name)_longitude"][:], ds["$(name)_latitude"][:], ds["$(name)_depth"][:]))
        return (; name, model, values=ds[sample][:], σ=ds["$(σname)_sample"][:], locations)
    end, names)
end, files.observations)
# Zonally averaged interior sources at iteration zero, the final iteration, and in the reference TMI.
modeled_names = Set(getproperty.(observation_specs, :model))
reference_names = NCDataset(ds -> Set(keys(ds)), files.reference)
source_names = NCDataset(ds -> begin
    [String(name) for name in keys(ds) if dimnames(ds[name]) == ("lon", "lat", "depth") && name in reference_names && !(name in modeled_names)]
end, files.initial)
if !isempty(source_names)
    fig, ax = plt.subplots(length(source_names); ncols=4, squeeze=false, figsize=(24, 4.5length(source_names)), sharex=true, sharey=true, constrained_layout=true)
    for (row, name) in enumerate(source_names)
        sources = map(file -> readsource(file, name, γ), (files.initial, files.final, files.reference))
        zonal_sources = map(source -> zonalaverage(Field(source.tracer, γ, source.name, source.longname, source.units)), sources)
        difference = zonal_sources[1] .- zonal_sources[2]
        limits = extrema(filter(isfinite, vcat(zonal_sources...)))
        difference_limit = max(maximum(abs, filter(isfinite, difference)), eps())
        titles = ("Iteration 0", "Iteration $final_iteration", reference, "Iteration 0 − iteration $final_iteration")
        maps = []
        for (column, (title, values)) in enumerate(zip(titles, (zonal_sources..., difference)))
            isdifference = column == 4
            levels = isdifference ? range(-difference_limit, difference_limit; length=21) : range(limits...; length=21)
            axis = ax[row - 1, column - 1]
            push!(maps, axis.contourf(γ.lat, -γ.depth, values'; levels, cmap=isdifference ? "RdBu_r" : nothing))
            axis.set(title="$(sources[1].longname) — $title", xlabel="Latitude", ylabel=column == 1 ? "Depth [m]" : "")
        end
        fig.colorbar(maps[3]; ax=[ax[row - 1, column] for column in 0:2], label=sources[1].units)
        fig.colorbar(maps[4]; ax=ax[row - 1, 3], label="$(sources[1].units) difference")
    end
    fig.savefig(plotfile("source_zonal_average_comparison"); dpi=200); plt.close(fig)
end
checkpoint_files = filter(file -> startswith(basename(file), "$(experiment)_iteration_") && endswith(file, ".nc"), readdir(checkpointdir; join=true))
progress = [NCDataset(ds -> begin
    misfits = map(function (observation)
        modeled = readfield(file, observation.model, γ; name=Symbol(observation.model))
        prediction = isnothing(observation.locations) ? modeled.tracer[γ.interior] : observe(modeled, observation.locations, γ)
        sqrt(mean(abs2, (prediction .- observation.values) ./ observation.σ))
    end, observation_specs)
    (file=file, iteration=Int(ds.attrib["iteration"]), objective=Float64(ds.attrib["objective"]), misfits)
end, file) for file in [checkpoint_files; files.final]]
progress = unique(x -> x.iteration, sort(progress; by=x -> x.iteration))
iterations = getproperty.(progress, :iteration)
fig, ax = plt.subplots(2; figsize=(10, 9), sharex=true)
ax[0].plot(iterations, getproperty.(progress, :objective); marker="o"); ax[0].set(ylabel="Objective function", yscale="log")
for (index, observation) in enumerate(observation_specs)
    ax[1].plot(iterations, [x.misfits[index] for x in progress]; marker="o", label=observation.name)
end
ax[1].axhline(1; color="0.5", linestyle="--"); ax[1].set(xlabel="Iteration", ylabel="Normalized RMS misfit", xscale="symlog", yscale="log"); ax[1].legend(ncols=3)
ax[1].set_xlim(left=0)
ax[0].grid(true; which="both", alpha=0.3); ax[1].grid(true; which="both", alpha=0.3)
fig.savefig(plotfile("cost_history"); dpi=200); plt.close(fig)
# Volume-filled evolution and reference comparison.
targets = round.(Int, range(0, last(iterations); length=9))
snapshots = map(target -> progress[argmin(abs.(iterations .- target))], targets)
filled = [TMI.volumefilled(NaN, lu(TMI.watermassmatrix(x.file)), γ) for x in snapshots]; reference_filled = TMI.volumefilled(NaN, lu(TMI.watermassmatrix(files.reference)), γ)
filled_maps = map(x -> 10.0 .^ x.tracer, [filled; reference_filled])
labels = [["Iteration $(x.iteration)" for x in snapshots]; reference]
limits = extrema(filter(isfinite, vcat(filled_maps...))); log_color = plt.matplotlib.colors.LogNorm(vmin=limits[1], vmax=limits[2])
function volumemap(ax, field, title)
    ax.set(title=title, xlabel="Longitude", ylabel="Latitude")
    return ax.pcolormesh(γ.lon, γ.lat, field'; norm=log_color, shading="auto")
end
fig, ax = plt.subplots(2; ncols=5, figsize=(20, 8), sharex=true, sharey=true, constrained_layout=true)
maps = [volumemap(a, field, label) for (a, field, label) in zip(ax.flat, filled_maps, labels)]
fig.colorbar(last(maps); ax=ax, label="m³/m²")
fig.savefig(plotfile("volume_filled_maps"); dpi=200); plt.close(fig)
fig, ax = plt.subplots()
for (snapshot, result, color) in zip(snapshots, filled, range(0, 1; length=length(filled)))
    ax.plot(sort(filter(isfinite, result.tracer)); color=plt.get_cmap("winter")(color), label="Iteration $(snapshot.iteration)")
end
ax.plot(sort(filter(isfinite, reference_filled.tracer)); color="black", label=reference)
ax.set(xlabel="Sorted surface grid cell", ylabel=first(filled).units); ax.legend()
fig.savefig(plotfile("volume_filled_distributions"); dpi=200); plt.close(fig)
fig, ax = plt.subplots(1; ncols=3, figsize=(19, 5.5), constrained_layout=true)
volumemap(ax[0], filled_maps[end - 1], "Experiment: $experiment")
volume_map = volumemap(ax[1], filled_maps[end], reference)
fig.colorbar(volume_map; ax=[ax[0], ax[1]], label="m³/m²")
ax[2].plot(sort(filter(isfinite, last(filled).tracer)); color="black", label="Experiment")
ax[2].plot(sort(filter(isfinite, reference_filled.tracer)); color="grey", label=reference)
ax[2].set(title="Volume-filled distributions", xlabel="Sorted surface grid cell", ylabel=first(filled).units); ax[2].legend()
fig.savefig(plotfile("volume_filled_reference_comparison"); dpi=200); plt.close(fig)

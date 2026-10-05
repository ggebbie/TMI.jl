#=
Plot an inversion saved by runinversion. Run from the repository with Julia 1.13,
with GeoPythonPlot in the global environment, which comes first in the load path:
    JULIA_LOAD_PATH="@v1.13:.:@stdlib" julia +1.13 --project=@v1.13 scripts/inversion_diagnostics.jl full_inversion
Figures: observed, first-guess, and final tracers; zonal-mean sources; cost and
misfit history; volume filled by each surface cell.
=#
ENV["JULIA_CONDAPKG_BACKEND"] = "Null"
ENV["JULIA_PYTHONCALL_EXE"] = get(ENV, "PYTHON", "python")
ENV["MPLBACKEND"] = "Agg"
using GeoPythonPlot, LinearAlgebra, NCDatasets, Statistics, TMI
const plt = GeoPythonPlot.pyplot

experiment = isempty(ARGS) ? "full_inversion" : only(ARGS)
reference = "modern_90x45x33_G14_v2"
checkpointdir = TMI.pkgdir("checkpoints", experiment)
plotdir = TMI.pkgplotsdir(experiment)
files = (observations = joinpath(checkpointdir, "$(experiment)_observations.nc"),
    initial = joinpath(checkpointdir, "$(experiment)_iteration_000000.nc"),
    final = TMI.pkgdatadir("$(experiment)_final.nc"),
    reference = TMI.pkgdatadir("TMI_$reference.nc"))
mkpath(plotdir)
println("Saving figures to $plotdir")
γ = Grid(files.reference)
final_iteration = NCDataset(ds -> Int(ds.attrib["iteration"]), files.final)
plotfile(name) = joinpath(plotdir, "$(experiment)_$name.png")
readtracer(file, name) = readfield(file, name, γ; name=Symbol(name))

# observed fields are NaN where a tracer is not observed, which readfield rejects
function readobserved(name)
    tracer, units, longname = TMI._read3d(files.observations, name)
    return Field(tracer, γ, Symbol(name), longname, units)
end

"""
    zonalmean(c)

Volume-weighted zonal mean of a `Field`, as a (lat, depth) matrix.
"""
function zonalmean(c::Field)
    weights = ifelse.(γ.wet .& isfinite.(c.tracer), cellvolume(γ).tracer, 0.0)
    values = ifelse.(weights .> 0, c.tracer, 0.0)
    return dropdims(sum(values .* weights; dims=1) ./ sum(weights; dims=1); dims=1)
end

"""
    comparisonrow!(fig, ax, x, y, states, difference, titles, units, xlabel, ylabel)

Contour three states with shared color levels, and a difference of two of them
in a fourth panel.
"""
function comparisonrow!(fig, ax, x, y, states, difference, titles, units, xlabel, ylabel)
    limits = extrema(filter(isfinite, vcat(states...)))
    difference_limit = max(maximum(abs, filter(isfinite, difference)), eps())
    maps = []
    for (column, values) in enumerate((states..., difference))
        isdifference = column == 4
        levels = isdifference ? range(-difference_limit, difference_limit; length=21) :
            range(limits...; length=21)
        axis = ax[column - 1]
        push!(maps, axis.contourf(x, y, values'; levels, cmap=isdifference ? "RdBu_r" : nothing))
        axis.set(title=titles[column], xlabel=xlabel, ylabel=column == 1 ? ylabel : "")
    end
    fig.colorbar(maps[3]; ax=[ax[column] for column in 0:2], label=units)
    fig.colorbar(maps[4]; ax=ax[3], label="$units difference")
end

# observed tracers: gridded values, or values and locations at sparse points
observations = NCDataset(files.observations) do ds
    names = filter(name -> haskey(ds[name].attrib, "support"), keys(ds))
    map(names) do name
        if ds[name].attrib["support"] == "gridded"
            values = ds[name][:, :, :][γ.interior]
            σ = ds["σ$name"][:, :, :][γ.interior]
            locations = nothing
        else
            values = ds["$(name)_sample"][:]
            σ = ds["σ$(name)_sample"][:]
            locations = collect(zip(ds["$(name)_longitude"][:], ds["$(name)_latitude"][:],
                ds["$(name)_depth"][:]))
        end
        (; name, values, σ, locations)
    end
end

# observed, first-guess, and final tracers at the surface and zonally averaged
surface = argmin(abs.(γ.depth))
for observation in observations
    fields = (readobserved(observation.name), readtracer(files.initial, observation.name),
        readtracer(files.final, observation.name))
    label = fields[2].longname
    titles = ("Observations", "Iteration 0", "Iteration $final_iteration", "Final − iteration 0")
    fig, ax = plt.subplots(2; ncols=4, figsize=(22, 9), sharex="row", sharey="row",
        squeeze=false, constrained_layout=true)
    surface_states = map(c -> c.tracer[:, :, surface], fields)
    comparisonrow!(fig, ax[0], γ.lon, γ.lat, surface_states, surface_states[3] .- surface_states[2],
        map(t -> "$label surface — $t", titles), fields[2].units, "Longitude", "Latitude")
    zonal_states = map(zonalmean, fields)
    comparisonrow!(fig, ax[1], γ.lat, -γ.depth, zonal_states, zonal_states[3] .- zonal_states[2],
        map(t -> "$label zonal mean — $t", titles), fields[2].units, "Latitude", "Depth [m]")
    fig.savefig(plotfile("$(observation.name)_quick_diagnostics"); dpi=200)
    plt.close(fig)
end

# zonal-mean interior sources: first guess, final, and reference TMI
modeled_names = Set(getproperty.(observations, :name))
reference_names = NCDataset(ds -> Set(keys(ds)), files.reference)
source_names = NCDataset(files.initial) do ds
    filter(name -> dimnames(ds[name]) == ("lon", "lat", "depth") &&
        name in reference_names && !(name in modeled_names), collect(keys(ds)))
end
if !isempty(source_names)
    fig, ax = plt.subplots(length(source_names); ncols=4, squeeze=false,
        figsize=(24, 4.5length(source_names)), sharex=true, sharey=true, constrained_layout=true)
    for (row, name) in enumerate(source_names)
        q = map(file -> readsource(file, name, γ), (files.initial, files.final, files.reference))
        zonal = map(qₖ -> zonalmean(Field(qₖ.tracer, γ, qₖ.name, qₖ.longname, qₖ.units)), q)
        titles = map(t -> "$(q[1].longname) — $t", ("Iteration 0", "Iteration $final_iteration",
            reference, "Iteration 0 − iteration $final_iteration"))
        comparisonrow!(fig, ax[row - 1], γ.lat, -γ.depth, zonal, zonal[1] .- zonal[2], titles,
            q[1].units, "Latitude", "Depth [m]")
    end
    fig.savefig(plotfile("source_zonal_average_comparison"); dpi=200)
    plt.close(fig)
end

# cost and normalized RMS misfit of each tracer at every saved iteration
checkpoint_files = filter(file -> startswith(basename(file), "$(experiment)_iteration_") &&
    endswith(file, ".nc"), readdir(checkpointdir; join=true))
progress = map([checkpoint_files; files.final]) do file
    misfits = map(observations) do observation
        c = readtracer(file, observation.name)
        ỹ = isnothing(observation.locations) ? c.tracer[γ.interior] :
            observe(c, observation.locations, γ)
        sqrt(mean(abs2, (ỹ .- observation.values) ./ observation.σ))
    end
    NCDataset(file) do ds
        (; file, iteration=Int(ds.attrib["iteration"]), J=Float64(ds.attrib["J"]),
            Jdata=Float64(ds.attrib["Jdata"]), Jcontrol=Float64(ds.attrib["Jcontrol"]), misfits)
    end
end
progress = unique(p -> p.iteration, sort(progress; by=p -> p.iteration))
iterations = getproperty.(progress, :iteration)
fig, ax = plt.subplots(2; figsize=(10, 9), sharex=true)
for term in (:J, :Jdata, :Jcontrol)
    ax[0].plot(iterations, getproperty.(progress, term); marker="o", label=String(term))
end
ax[0].set(ylabel="Cost function", yscale="log")
ax[0].legend()
for (i, observation) in enumerate(observations)
    ax[1].plot(iterations, [p.misfits[i] for p in progress]; marker="o", label=observation.name)
end
ax[1].axhline(1; color="0.5", linestyle="--")
ax[1].set(xlabel="Iteration", ylabel="Normalized RMS misfit", xscale="symlog", yscale="log")
ax[1].legend(ncols=3)
ax[1].set_xlim(left=0)
ax[0].grid(true; which="both", alpha=0.3)
ax[1].grid(true; which="both", alpha=0.3)
fig.savefig(plotfile("cost_history"); dpi=200)
plt.close(fig)

# volume filled by each surface cell, through the iterations and in the reference TMI
targets = round.(Int, range(0, last(iterations); length=9))
snapshots = map(target -> progress[argmin(abs.(iterations .- target))], targets)
filled = [TMI.volumefilled(NaN, lu(TMI.watermassmatrix(p.file)), γ) for p in snapshots]
reference_filled = TMI.volumefilled(NaN, lu(TMI.watermassmatrix(files.reference)), γ)
filled_maps = map(v -> 10.0 .^ v.tracer, [filled; reference_filled])
labels = [["Iteration $(p.iteration)" for p in snapshots]; reference]
limits = extrema(filter(isfinite, vcat(filled_maps...)))
log_color = plt.matplotlib.colors.LogNorm(vmin=limits[1], vmax=limits[2])
function volumemap(axis, values, title)
    axis.set(title=title, xlabel="Longitude", ylabel="Latitude")
    return axis.pcolormesh(γ.lon, γ.lat, values'; norm=log_color, shading="auto")
end
fig, ax = plt.subplots(2; ncols=5, figsize=(20, 8), sharex=true, sharey=true, constrained_layout=true)
maps = [volumemap(axis, values, label) for (axis, values, label) in zip(ax.flat, filled_maps, labels)]
fig.colorbar(last(maps); ax=ax, label="m³/m²")
fig.savefig(plotfile("volume_filled_maps"); dpi=200)
plt.close(fig)

fig, ax = plt.subplots()
for (snapshot, v, color) in zip(snapshots, filled, range(0, 1; length=length(filled)))
    ax.plot(sort(filter(isfinite, v.tracer)); color=plt.get_cmap("winter")(color),
        label="Iteration $(snapshot.iteration)")
end
ax.plot(sort(filter(isfinite, reference_filled.tracer)); color="black", label=reference)
ax.set(xlabel="Sorted surface grid cell", ylabel=first(filled).units)
ax.legend()
fig.savefig(plotfile("volume_filled_distributions"); dpi=200)
plt.close(fig)

fig, ax = plt.subplots(1; ncols=3, figsize=(19, 5.5), constrained_layout=true)
volumemap(ax[0], filled_maps[end - 1], "Experiment: $experiment")
volume_map = volumemap(ax[1], filled_maps[end], reference)
fig.colorbar(volume_map; ax=[ax[0], ax[1]], label="m³/m²")
ax[2].plot(sort(filter(isfinite, last(filled).tracer)); color="black", label="Experiment")
ax[2].plot(sort(filter(isfinite, reference_filled.tracer)); color="grey", label=reference)
ax[2].set(title="Volume-filled distributions", xlabel="Sorted surface grid cell",
    ylabel=first(filled).units)
ax[2].legend()
fig.savefig(plotfile("volume_filled_reference_comparison"); dpi=200)
plt.close(fig)


import Pkg; Pkg.activate("./scripts")
Pkg.instantiate()
using Revise
using PlutoUI
using TMI
using Test
using PythonCall
using CondaPkg
using PythonPlot
const cartopy = pyimport("cartopy")
const matplotlib = pyimport("matplotlib")

ccrs = cartopy.crs

import Pkg;
using Interpolations
using Statistics
using LinearAlgebra
plotsdir(x) = "/Users/anthonymeza/Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology/Documents/GitHub/TMI.jl/plots/" * x
include("../src/sample_observations.jl")
# take synthetic, noisy observations

wetsurfacelocation(γ)

A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);
Nd = length(γ.depth)
N = 20
function synthetic_benthic_and_surface_observation(TMIversion; N = 25, σ = 0)
	y, W⁻, θtrue, ytrue, locs, wis = synthetic_column_observations(TMIversion,"θ",γ,N, σ);

	bottom_indices = [findfirst(isnan, y[:, i]) for i in 1:N];
	bottom_indices[(!isnothing).(bottom_indices)] .-= 1;
	bottom_indices[isnothing.(bottom_indices)] .= Nd;

	benthic_observations = [y[bi, i] for (i, bi) in enumerate(bottom_indices)];
	surface_observations = [y[1, i] for (i, bi) in enumerate(bottom_indices)];

	return (benthic_observations, surface_observations, locs)
end

function synthetic_benthic_and_surface_observation(TMIversion, locs ; σ = 0)
	N = length(locs)
	y, W⁻, θtrue, ytrue, locs, wis = synthetic_column_observations(TMIversion,"θ",γ, locs, σ);

	bottom_indices = [findfirst(isnan, y[:, i]) for i in 1:N];
	bottom_indices[(!isnothing).(bottom_indices)] .-= 1;
	bottom_indices[isnothing.(bottom_indices)] .= Nd;

	benthic_observations = [y[bi, i] for (i, bi) in enumerate(bottom_indices)];
	surface_observations = [y[1, i] for (i, bi) in enumerate(bottom_indices)];

	return (benthic_observations, surface_observations, locs)
end

niter = 50
benthic_difference = zeros(niter)
surface_difference = zeros(niter)
average_latitude = zeros(niter)

for i in 1:niter
	TMIversion = "modern_90x45x33_G14_v2";

	modern_observations = Dict()
	benthic_observations, surface_observations, locs = synthetic_benthic_and_surface_observation(TMIversion; N = 100, σ = 0);

	modern_observations["benthic"] = benthic_observations
	modern_observations["surface"] = surface_observations

	TMIversion = "LGM_90x45x33_G14";
	LGM_observations = Dict()
	benthic_observations, surface_observations, locs = synthetic_benthic_and_surface_observation(TMIversion, locs; σ = 0);
	LGM_observations["benthic"] = benthic_observations
	LGM_observations["surface"] = surface_observations

	benthic_difference[i] = mean(modern_observations["benthic"] .- LGM_observations["benthic"])
	surface_difference[i] = mean(modern_observations["surface"] .- LGM_observations["surface"])

	average_latitude[i] = mean([l[2] for l in locs])
end

surface_average(x, γ) = sum(x.tracer[wet(x)] .* cellarea(γ).tracer[wet(x)]) / sum(cellarea(γ).tracer[wet(x)])

TMIversion = "modern_90x45x33_G14_v2";
TMIfile = TMI.pkgdatadir("TMI_"*TMIversion*".nc")
θ̄_mod = readfield(TMIfile, "θ", γ); #true mean
θ̄_mod_MOT = mean(θ̄_mod); #true mean
θ̄_mod_SST = surface_average(getsurfaceboundary(θ̄_mod), γ);

TMIversion = "LGM_90x45x33_G14";
TMIfile = TMI.pkgdatadir("TMI_"*TMIversion*".nc")
θ̄_LGM = readfield(TMIfile, "θ", γ); #true mean
θ̄_LGM_MOT = mean(θ̄_LGM); #true mean
θ̄_LGM_SST = surface_average(getsurfaceboundary(θ̄_LGM), γ);

fig, ax = subplots(figsize = (6, 5))
ax.plot(collect(0:5), collect(0:5), linestyle = "--", c = "k")
ax.scatter(surface_difference, benthic_difference, alpha = 0.2, edgecolor = "grey")
ax.scatter(θ̄_mod_SST - θ̄_LGM_SST, θ̄_mod_MOT - θ̄_LGM_MOT, s= 200, edgecolor = "k")
ax.set_xlim(0, 5)
ax.set_ylim(0, 5)
ax.set_xlabel("ΔSST")
ax.set_ylabel("ΔMOT")
fig

fig, ax = subplots(figsize = (6, 5))
ax.plot(collect(0:5), collect(0:5), alpha = 0.5, linestyle = "--", c = "k")
cmp = ax.scatter(surface_difference, benthic_difference, alpha = 0.9, edgecolor = "grey", 
		   c = average_latitude, cmap = "bwr", vmin = -65, vmax = 65)
fig.colorbar(cmp, label = "Average Latitude")
ax.scatter(θ̄_mod_SST - θ̄_LGM_SST, θ̄_mod_MOT - θ̄_LGM_MOT, s= 200, edgecolor = "k")
ax.set_xlim(0, 5)
ax.set_ylim(0, 5)
ax.set_xlabel("ΔSST", fontweight = "bold")
ax.set_ylabel("ΔMOT", fontweight = "bold")
fig


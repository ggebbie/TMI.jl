
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

versionlist()
TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);


volume = volumefilled(TMIversion,Alu,γ)
volume = (10 .^ volume.tracer) .* cellarea(γ).tracer
volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","m³")

figl, axl = subplots(figsize = (10, 10), 
						subplot_kw=Dict("projection"=> ccrs.PlateCarree()))

# Normalize this manually (no need for LogNorm now)
norm = matplotlib.colors.LogNorm(vmin = 2e11, vmax=2e17)
axl.coastlines(color = "#949494", alpha = 0.3)
cml = axl.pcolormesh(γ.lon, γ.lat, volume.tracer', cmap = "bwr", transform = ccrs.PlateCarree(),norm=norm)
cb = figl.colorbar(cml, orientation = "horizontal", extend = "both")
figl

volume = volumefilled(TMIversion,Alu,γ)
volume = (10 .^ volume.tracer) .* cellarea(γ).tracer

volume = 100 .* volume / sum(cellvolume(γ).tracer)
volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","m³")


figl, axl = subplots(figsize = (10, 10), 
						subplot_kw=Dict("projection"=> ccrs.PlateCarree()))
norm = matplotlib.colors.LogNorm(vmin = 1e-5, vmax=100)

# Normalize this manually (no need for LogNorm now)
axl.coastlines(color = "#949494", alpha = 0.3)
cml = axl.pcolormesh(γ.lon, γ.lat, volume.tracer', cmap = "bwr", 
transform = ccrs.PlateCarree(), norm = norm)
cb = figl.colorbar(cml, orientation = "horizontal", extend = "both")
figl



TMIversion = "modern_90x45x33_G14_v2";
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion); println(sum(γ.wet))
volume = volumefilled(TMIversion,Alu,γ)
volume = (10 .^ volume.tracer) .* cellarea(γ).tracer
volume = volume / sum(cellvolume(γ).tracer)
volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","unknown")
volume_mod = 1 * volume 
TMIfile = TMI.pkgdatadir("TMI_"*TMIversion*".nc")
θ̄_mod = readfield(TMIfile, "θ", γ); #true mean
θ̄_mod_MOT = mean(θ̄_mod); #true mean
θ̄_mod_SST = getsurfaceboundary(θ̄_mod) * volume; θ̄_mod_SST = sum(θ̄_mod_SST.tracer[γ.wet[:, :, 1]])

TMIversion = "LGM_90x45x33_G14";
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion); println(sum(γ.wet))
volume = volumefilled(TMIversion,Alu,γ)
volume = (10 .^ volume.tracer) .* cellarea(γ).tracer
volume = volume / sum(cellvolume(γ).tracer)
volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","unknown")
volume_LGM = 1 * volume
TMIfile = TMI.pkgdatadir("TMI_"*TMIversion*".nc")
θ̄_LGM = readfield(TMIfile, "θ", γ); #true mean
θ̄_LGM_MOT = mean(θ̄_LGM); #true mean
θ̄_LGM_SST = getsurfaceboundary(θ̄_LGM) * volume; θ̄_LGM_SST = sum(θ̄_LGM_SST.tracer[γ.wet[:, :, 1]])

θ̄_mod_SST - θ̄_LGM_SST



ΔθW = θ̄_mod_SST - θ̄_LGM_SST 

Δθ = getsurfaceboundary(θ̄_mod - θ̄_LGM); 

Δv_norm = sqrt(sum(volume_LGM.tracer[γ.wet[:, :, 1]].^2))
Δθ_norm = sqrt(sum(Δθ.tracer[γ.wet[:, :, 1]].^2))

Sc = ΔθW / (Δv_norm * Δθ_norm) #cosine similarity, near zero means orthogonal
rad2deg(acos(Sc))
#orthogonality means that these two pieces carry unique information ? 

mean(Δθ)
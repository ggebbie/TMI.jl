
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

TMIversion = "modern_90x45x33_G14_v2";
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion); println(sum(γ.wet))
volume = volumefilled(TMIversion,Alu,γ)
volume = (10 .^ volume.tracer) .* cellarea(γ).tracer
volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","unknown")
volume_mod = 1 * volume 
TMIfile = TMI.pkgdatadir("TMI_"*TMIversion*".nc")
θ̄_mod = readfield(TMIfile, "θ", γ); #true mean
θ̄_mod_MOT = mean(θ̄_mod); #true mean
θ̄_mod_SST = getsurfaceboundary(θ̄_mod)

TMIversion = "LGM_90x45x33_G14";
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion); println(sum(γ.wet))
volume = volumefilled(TMIversion,Alu,γ)
volume = (10 .^ volume.tracer) .* cellarea(γ).tracer
volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","unknown")
volume_LGM = 1 * volume
TMIfile = TMI.pkgdatadir("TMI_"*TMIversion*".nc")
θ̄_LGM = readfield(TMIfile, "θ", γ); #true mean
θ̄_LGM_MOT = mean(θ̄_LGM); #true mean
θ̄_LGM_SST = getsurfaceboundary(θ̄_LGM)

figl, ax = subplots(2, 3, figsize = (10, 10), 
						subplot_kw=Dict("projection"=> ccrs.Robinson(central_longitude = 180)))

labels = ["LGM Volume Contribution", "Modern Volume Contribution", "Difference\n(LGM minus Modern)"] 
volumes = [volume_LGM.tracer, volume_mod.tracer, -(volume_mod.tracer .- volume_LGM.tracer)]
vnorms = [matplotlib.colors.LogNorm(vmin = 2e11, vmax=2e16), 
matplotlib.colors.LogNorm(vmin = 2e11, vmax=2e16), matplotlib.colors.SymLogNorm(linthresh =1e14, vmin = -2e15, vmax=2e15)]

vnorm = vnorms[1]
cmaps = ["Spectral_r", "Spectral_r", "bwr"]

for (i, axl) in enumerate(ax[0, [0, 1, 2]])
    axl.set_title(labels[i])
    axl.coastlines(color = "#949494", alpha = 0.3)
    cml = axl.pcolormesh(γ.lon, γ.lat, volumes[i]', cmap = cmaps[i], transform = ccrs.PlateCarree(),norm=vnorms[i])
    cb = figl.colorbar(cml, orientation = "horizontal", 
    extend = "both", label = L"[m^3]")
end
figl

labels = ["LGM SST", "Modern SST", "Difference\n(LGM minus Modern)"] 

temperatures = [θ̄_LGM_SST.tracer, θ̄_mod_SST.tracer, -(θ̄_mod_SST.tracer .- θ̄_LGM_SST.tracer)]
vnorms = [matplotlib.colors.Normalize(vmin = -2, vmax=30), 
matplotlib.colors.Normalize(vmin = -2, vmax=30), matplotlib.colors.Normalize(vmin = -10, vmax=10)]
cmaps = ["magma", "magma", "bwr"]
for (i, axl) in enumerate(ax[1, [0, 1, 2]])
    axl.set_title(labels[i])
    axl.coastlines(color = "#949494", alpha = 0.3)
    cml = axl.pcolormesh(γ.lon, γ.lat, temperatures[i]', 
    cmap = cmaps[i], 
    transform = ccrs.PlateCarree(),norm=vnorms[i])
    cb = figl.colorbar(cml, orientation = "horizontal", extend = "both", label = "[deg C]")
end
figl.tight_layout()
figl.savefig(plotsdir("SST_Volume_Differences.png"))

# volume = volumefilled(TMIversion,Alu,γ)
# volume = (10 .^ volume.tracer) .* cellarea(γ).tracer

# volume = 100 .* volume / sum(cellvolume(γ).tracer)
# volume = BoundaryCondition(volume,(γ.lon,γ.lat),γ.depth[1],3,1,γ.wet[:,:,1],:V,"volume filled by surface gridcell","m³")


# figl, axl = subplots(figsize = (10, 10), 
# 						subplot_kw=Dict("projection"=> ccrs.PlateCarree()))
# norm = matplotlib.colors.LogNorm(vmin = 1e-5, vmax=100)

# # Normalize this manually (no need for LogNorm now)
# axl.coastlines(color = "#949494", alpha = 0.3)
# cml = axl.pcolormesh(γ.lon, γ.lat, volume.tracer', cmap = "bwr", 
# transform = ccrs.PlateCarree(), norm = norm)
# cb = figl.colorbar(cml, orientation = "horizontal", extend = "both")
# figl
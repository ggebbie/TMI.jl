
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

ΔθW = θ̄_mod_MOT - θ̄_LGM_MOT 

Δθ = getsurfaceboundary(θ̄_mod - θ̄_LGM); 
WΔθ = Δθ * volume_LGM; WΔθ = sum(WΔθ.tracer[γ.wet[:, :, 1]])

ΔW = (volume_mod - volume) 
θΔW = ΔW * getsurfaceboundary(θ̄_LGM); θΔW = sum(θΔW.tracer[γ.wet[:, :, 1]])

ΔWΔθ = Δθ * ΔW; ΔWΔθ = sum(ΔWΔθ.tracer[γ.wet[:, :, 1]])

fig, ax = subplots(figsize = (5, 7.5))
groups = [L"v_{LGM} \Delta \theta", L"\theta_{LGM} \Delta v", L"\Delta v\Delta \theta"]

bar_width = 0.1
b1= ax.bar(1, WΔθ, color = "r", label = [groups[1]], width=bar_width, alpha = 0.7)
b2 = ax.bar(1, θΔW, bottom=WΔθ, color = "b", label = [groups[2]], width=bar_width, alpha = 0.7)
b3 = ax.bar(1, ΔWΔθ, bottom= WΔθ + θΔW, color = "g", label = [groups[3]], width=bar_width, alpha = 0.7)
# ax.legend()
fig
# Add labels to each section of the bar
for bars in [b1, b2, b3]
    for bar in bars
        yval = bar.get_height()
        yval_r = round(pyconvert(Float32, yval), digits = 2)
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_y() + yval/2, string(yval_r) * " deg C", 
                ha="center", va="center")
    end
end
fig

# Calculate the total height of the stacked bar
total_height = 1 * round(ΔθW, digits = 2)

# Add horizontal curly bracket along the side of the bar using ax.transAxes
ax.annotate("", xy=(.7, 0.1), xytext=(0.7, 0.9),
            arrowprops=Dict("arrowstyle"=>"-", "lw"=>2),
            xycoords="axes fraction", textcoords="axes fraction")

# Add total label next to the bar with rotation
ax.text(0.75, 0.5, "Total MOT Change: " * string(total_height) * " deg C", va="center", rotation=90, transform=ax.transAxes)

# Customize plot
ax.set_xticks([1])
ax.set_xticklabels(["MOT Change"])
ax.set_xlim(0.8, 1.2)  # Adjust x-axis limits to fit the skinny bar and bracket
ax.set_ylabel("[deg C]")
ax.set_title("Difference in MOT between Modern Day and LGM")

# Show legend
ax.legend()

# Show plot
fig.savefig(plotsdir("MOT_Change_decomposition.png"), dpi = 200)
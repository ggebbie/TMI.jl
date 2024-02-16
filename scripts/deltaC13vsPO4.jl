#=
% Example: Compare average δ¹³C and PO₄ values for different ocean regions %
=#
# import Pkg 
# Pkg.activate("TMI.jl")#do this if using script interactively

using Revise
using TMI, PyPlot, PyCall, Statistics, DrWatson, PyCall, Colors

@pyimport matplotlib.colors as mcolors
@pyimport matplotlib.colors as plt_colors
mkpath(plotsdir())

TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
list = ("NATL","NPAC","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

# contains the percentage of volume sourced from each 
# surface ocean 
depth = getVar(TMIversion, "depth")

#initialize a grid containing the water mass sources
sources = identifywatermass(TMIversion,Alu,γ, list)

lims = 1:1.0:length(list)

cm = get_cmap(:tab20)
colorrange = (0:(length(list)-1)) ./ (length(list))

XC = γ.lon' .* ones(length(γ.lat))
YC = ones(length(γ.lon))' .*  γ.lat

level = 20;
sources_in_layer = sources[:, :, level] 
source_values = sort(filter(!isnan, unique(sources_in_layer)))
nsource = length(source_values)

#making a color bar
cmap = cm(colorrange[Int.(source_values)]); 
my_cmap = plt_colors.ListedColormap(cmap)
bounds = source_values # define the bins and normalize
norm = plt_colors.BoundaryNorm(bounds, my_cmap.N)

#plot the distribution of source waters within a layer
fig, ax = subplots(1,1)
good_pts = (!isnan).(sources_in_layer) # filter for nonnan values
[ax.scatter(XC'[good_pts][1], YC'[good_pts][1], c = val, cmap = my_cmap, label = list[val], norm=norm, zorder = 0) for val in source_values] #plot the labels
fig.legend(ncols = 3)
ax.scatter(XC'[good_pts][:],YC'[good_pts][:], c = sources_in_layer[good_pts][:], cmap = my_cmap, norm=norm) #plot all points
fig.savefig(plotsdir("RegionalMasks.png"))
fig

#Plot PO vs C13 for various regions, shade the dots by their source region
fig, ax = subplots(1,1)
basin_list = ["","SUBANTIND", "NPAC",
            "SUBANTPAC","NATL","SUBANTATL"]
δ¹³C = getVar(TMIversion,"δ¹³C")
PO₄ = getVar(TMIversion,"PO₄")
fig, axs = plt.subplots(2, 3, figsize = (13, 8), sharex = true, sharey = true)
nb = length(basin_list)
[a.set_xlabel("PO₄") for a in axs]
[a.set_ylabel("δ¹³C") for a in axs]

for i in 2:nb
    ax = axs[i] 
    basin_mask = getRegion(TMIversion,"d_" * basin_list[i])

    ax.set_title(basin_list[i])
    masked_sources = sources[:, :, level] .* basin_mask
    masked_sources[iszero.(masked_sources)] .= NaN

    good_pts = (!isnan).(mask_water_mass)
    println("good points are: ", size(good_pts), "(should be geq to og + addtl pts)")
    masked_sources = masked_sources[good_pts]
    mask_δ¹³C = (δ¹³C[:, :, level] .* basin_mask)[good_pts]
    mask_PO₄ = (PO₄[:, :, level] .* basin_mask)[good_pts]
    ax.scatter(mask_PO₄[:],mask_δ¹³C[:], c = masked_sources[:], cmap = my_cmap, norm=norm)
end

depth_range = string(depth[level])
fig.suptitle("PO₄vsδ¹³C for various regions at "  * depth_range * "m")
fig.tight_layout()
fig.savefig("PO₄vsδ¹³C_" * depth_range * ".png")
fig
#=
Create maps of origin of specified core locations and take difference. 
Following ex5
Uses data from Common Era proposal - brynn currently has these csvs 
=#


using TMI
using Revise
using PyPlot, PyCall, CSV, DataFrames, NaNMath

pygui(true) #needed for Atom, not sure what it will do in other places
filepath = "/home/brynn/Documents/JP/rdata/coreloc1.csv"
x = CSV.read(filepath,DataFrame,header = 1)
lats = x[:,2]
lons = x[:,3]
#convert to degrees East
lons = 360 .- lons
depths = x[:,4]
num = x[:,1]

ccrs = pyimport("cartopy.crs")

#scatter plot core locations
figure()
ax = subplot(projection = ccrs.PlateCarree())
ax.coastlines()
scatter(lons, lats,transform=ccrs.PlateCarree())
for i in 1:length(lons)
    text(x = lons[i], y = lats[i], s = string(i), size = "medium", transform = ccrs.PlateCarree())
end
ax.gridlines()
ax.set_extent([290, 350, 50, 70])
#savefig("images/coreloc1_abbrev.png", dpi = 300)

#from ex5.surfaceorigin.jl
TMIversion = "modern_90x45x33_GH10_GH12"
TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
#for every position, plot surface origin
for i in 1:length(lons)
    xlon = lons[i]
    xlat = lats[i]
    xdepth = depths[i]
    loc = (xlon,xlat,xdepth)
    origin = surfaceorigin(loc, Alu, γ)

    # units: effective thickness in log10(meters)
    surfacedensity(γ.lon,γ.lat,log10.(origin'), xlon, xlat, "images/" * string(num[i]) *".png")
end


#helper function: calculates surface origin of given indexed core location 
function returnorigin(i)
    xlon = lons[i]
    xlat = lats[i]
    xdepth = depths[i]
    loc = (xlon,xlat,xdepth)
    origin = surfaceorigin(loc, Alu, γ)
    return origin
end


using Statistics
g1_ind = 1:7
g2_ind = 8:length(num)

g1 = zeros((length(g1_ind), length(γ.lon), length(γ.lat)))
[g1[i, :, :,] = returnorigin(g1_ind[i]) for i ∈ 1:length(g1_ind)]
origin1 = mean(g1, dims = 1)[1,:,:]

g2 = zeros((length(g2_ind), length(γ.lon), length(γ.lat)))
[g2[i, :, :,] = returnorigin(g2_ind[i]) for i ∈ 1:length(g2_ind)]
origin2 = mean(g2, dims = 1)[1,:,:]

#Make surface density maps and save 
surfacedensity(γ.lon,γ.lat,log10.(origin1'), 0, 0, "images/mean1.png")
surfacedensity(γ.lon,γ.lat,log10.(origin2'), 0, 0, "images/mean2.png")
#take difference between the two - some nonsense with TMI.jl to make this work 
difference = origin1 - origin2
difference = log10.(abs.(difference))
surfacedensity(γ.lon,γ.lat,difference', 0, 0, "images/diff.png")

difference = origin1 - origin2


figure()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.PlateCarree())
ax.coastlines()

pc = pcolormesh(γ.lon,γ.lat,difference', cmap = "coolwarm", vmin = -0.04, vmax = 0.04)
cb = colorbar(pc, fraction = 0.018)

cb.set_label("Effective Thickness [m]")

#add gridlines
gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
gl.top_labels = false
gl.right_labels = false

ax.set_extent([290, 20, 50, 80])

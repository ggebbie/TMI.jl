using TMI

using Revise
using PyPlot, PyCall, CSV, DataFrames

pygui(true) #needed for Atom, not sure what it will do in other places

filepath = "/home/brynn/Documents/JP/rdata/coreloc1_abbrev.csv"
x = CSV.read(filepath,DataFrame,header = 1)
lats = x[:,2]
lons = x[:,3]
#convert to degrees East
lons = 360 .- lons
depths = x[:,4]

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
ax.set_extent([330, 350, 50, 70])
savefig("images/coreloc1_abbrev.png", dpi = 300)

#from ex5.surfaceorigin.jl
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
#for every position, plot surface origin
for i in 1:length(lons)
    xlon = lons[i]
    xlat = lats[i]
    xdepth = depths[i]
    loc = (xlon,xlat,xdepth)
    origin = surfaceorigin(loc, Alu, γ)

    # units: effective thickness in log10(meters)
    surfacedensity(γ.lon,γ.lat,log10.(origin'), xlon, xlat, "images/ce_origin/" * string(i) *".png")
end

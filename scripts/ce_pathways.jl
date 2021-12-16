#compare pathways of (1,5) and (3,4)

using TMI

using Revise
using PyPlot, PyCall, CSV, DataFrames, Statistics

pygui(true) #needed for Atom, not sure what it will do in other places

filepath = "/home/brynn/Documents/JP/rdata/coreloc1_abbrev.csv"
x = CSV.read(filepath,DataFrame,header = 1)
lats = x[:,2]
lons = x[:,3]
#convert to degrees East
lons = 360 .- lons
depths = x[:,4]

lat1 = mean([lats[3], lats[4]])
lat2 = mean([lats[1], lats[5]])
lon1 = mean([lons[3], lons[4]])
lon2 = mean([lons[1], lons[5]])

ccrs = pyimport("cartopy.crs")

#scatter plot mean locations
figure()
ax = subplot(projection = ccrs.PlateCarree())
ax.coastlines()
scatter([lon1,lon2],[lat1,lat2],transform=ccrs.PlateCarree())
text(lon1, lat1, s = "1: " *string([lat1, lon1]), transform = ccrs.PlateCarree())
text(lon2, lat2, s = "2: "*string([lat2, lon2]), transform = ccrs.PlateCarree())
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
          linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
ax.set_extent([330, 350, 50, 70])
savefig("images/coreloc1_pathways.png", dpi = 300)

#do ex1.trackpathways, each location has a 1x1 degree box
#- define the surface patch by the bounding latitude and longitude.
# latbox1 = [60,62]
# lonbox1 = [335,337]
# latbox2 = [61,63]
# lonbox2 = [339,341]
latbox1= [58,62]
lonbox1=[334,338]
latbox2 = [58,62]
lonbox2=[338,342]

fname1 = ["images/ce_bbox1.png", "images/ce_dp1.png"]
fname2 = ["images/ce_bbox2.png", "images/ce_dp2.png"]
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

function track(latbox, lonbox, fnames)

    # do numerical analysis
    c = trackpathways(Alu,latbox,lonbox,γ)

    # do plotting (could be a function)
    plotextent(latbox, lonbox, fnames[1])

    # plot a section at the center of the longitudinal extent
    lon_section = lonbox[1]
    csection = section(c,lon_section,γ)
    println("shape = "*string(size(c)))
    lims = 0:5:100

    # make a plot of dye in the ocean
    dyeplot(γ.lat,-γ.depth[33:-1:1],100 * csection[:,33:-1:1]', lims, fnames[2])
end

track(latbox1, lonbox1, fname1)
track(latbox2, lonbox2, fname2)

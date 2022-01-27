#=
Create maps of origin of specified core locations and take difference. 
Following ex5
Uses data from Common Era proposal - brynn currently has these csvs 
=#
using TMI
using Revise
using PyPlot, PyCall, CSV, DataFrames, NaNMath, Statistics 
ccrs = pyimport("cartopy.crs")

#some helper functions
function surfacedensity(lon, lat, values, t) #could go in src 

    #init GeoAxes
    fig = figure()
    ax = fig.add_subplot(projection = ccrs.PlateCarree())
    ax.coastlines() #show coastlines
    
    cf = contourf(lon, lat, values, cmap = "plasma", levels = 50, vmax = 0, vmin = -10)
    cb = colorbar(cf,fraction = 0.018)
    cb.set_label("Effective Thickness [log10(m)]")

    #add gridlines
    gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
    gl.top_labels = false
    gl.right_labels = false
    title(t)

    ax.set_extent([290, 20, 50, 80])
end

#calculates surface origin of given indexed core location 
function returnorigin(i)
    xlon = lons[i]
    xlat = lats[i]
    xdepth = depths[i]
    loc = (xlon,xlat,xdepth)
    origin = surfaceorigin(loc, Alu, γ)
    return origin
end

#load in core location data - from CSV 
filepath = "/home/brynn/Documents/JP/rdata/coreloc1.csv"
x = CSV.read(filepath,DataFrame,header = 1)
lats = x[:,2]
lons = x[:,3]
#convert to degrees East
lons = 360 .- lons
depths = x[:,4]
num = x[:,1]

#TMIversion = "modern_90x45x33_GH10_GH12"
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
    surfacedensity(γ.lon,γ.lat,log10.(origin'), string(num[i]))
end

#we know that there are two distinct groups fom mse.jl 
g1_ind = 1:7
g2_ind = 8:length(num)

g1 = zeros((length(g1_ind), length(γ.lon), length(γ.lat)))
[g1[i, :, :,] = returnorigin(g1_ind[i]) for i ∈ 1:length(g1_ind)]
origin1 = Statistics.mean(g1, dims = 1)[1,:,:]

g2 = zeros((length(g2_ind), length(γ.lon), length(γ.lat)))
[g2[i, :, :,] = returnorigin(g2_ind[i]) for i ∈ 1:length(g2_ind)]
origin2 = mean(g2, dims = 1)[1,:,:]

#Make surface density maps and save 
surfacedensity(γ.lon,γ.lat,log10.(origin1'), "Mean of " *string(num[g1_ind]))
surfacedensity(γ.lon,γ.lat,log10.(origin2'), "Mean of " *string(num[g2_ind]))
#take difference between the two
difference = origin1 - origin2
difference = log10.(abs.(difference))
surfacedensity(γ.lon,γ.lat,difference', "ABS(LOG(DIFF))")

#now show difference without abs(log(
#something funky with cartopy here? can't use contourf
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

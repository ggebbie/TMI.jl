using PyPlot, CSV, DataFrames, PyCall
using TMI
filepath = "/home/brynn/Documents/JP/rdata/coreloc1.csv"
x = CSV.read(filepath,DataFrame,header = 1)
lats = x[:,2]
lons = x[:,3]
#convert to degrees East
lons = 360 .- lons
depths = x[:,4]
num = x[:, 1]


TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#land = sum(1 .- γ.depth .* γ.wet, dims = 3)[:,:,1]

#using PyCall
#ccrs = pyimport("cartopy.crs")
#ax = subplot(projection = ccrs.PlateCarree())
#lonlims = 150:180
#latlims = 70:80
#contourf(γ.lon[lonlims],γ.lat[latlims],land[lonlims,latlims]', cmap = "terrain",transform = ccrs.PlateCarree())
#scatter(lons,lats, s = 20, c = "red",transform = ccrs.PlateCarree())
#ax.coastlines()

figure()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection = ccrs.PlateCarree())

terrain = zeros((length(γ.lon), length(γ.lat)))
for i in 1:length(γ.lon)
    for j in 1:length(γ.lat)
        val = NaNMath.maximum(γ.wet[i,j,:] .* γ.depth)
        val = val < 1 ? NaN : val
        terrain[i,j] = val
    end
end


cf = pcolormesh(γ.lon, γ.lat, terrain', cmap = "Greens")
colorbar(cf)

s = scatter(lons,lats, s = 10, c = .- depths, cmap = "plasma", transform = ccrs.PlateCarree())
colorbar(s)

for i in 1:length(num)
    ax.text(x = lons[i], y = lats[i], s = string(num[i]),transform = ccrs.PlateCarree())
end


ax.coastlines()

gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=true,
                  linewidth=0.5, color="gray", alpha=0.5, linestyle="--")
gl.right_labels = false
gl.top_labels = false
#lat lon labels
ax.text(0.5, -0.12, "longitude", va="bottom", ha="center",
        rotation="horizontal", rotation_mode="anchor",
        transform=ax.transAxes)
ax.text(-0.12, 0.55, "latitude", va="bottom", ha="center",
        rotation="vertical", rotation_mode="anchor",
        transform=ax.transAxes)
ax.set_extent((335,341,59,66))


figure()
ax = gca()
scatter(lats,.-depths)
xlabel("latitude [°N]")
ylabel("depth [m]")
hlines(.- γ.depth, xmin = 58, xmax = 64, color = "green")
xlim((59,64))
ylim((-2500, -800))
for i in 1:length(num)
    ax.text(x = lats[i], y = .-  depths[i], s = string(num[i]))
end

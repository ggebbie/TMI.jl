using TMI, Revise, CSV, DataFrames, NCDatasets, PyPlot, Statistics
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

file = "/home/brynn/Documents/JP/rdata/ces_output.nc"

var = "theta"
data = NCDataset(file)[var][:]
age = NCDataset(file)["time"][:]

#get observations at each closest grid cell
ind = []
[push!(ind, (lons[i],lats[i],depths[i])) for i ∈ 1:length(lons)]
values = location_obs(data[:,:,:,:], ind, γ)
sind = []
[push!(sind, (lons[i], lats[i],0)) for i ∈ 1:length(lons)]
surface_values = location_obs(data, sind, γ)
surface_values = mean(surface_values, dims = 1)[1,:]

figure()
for i in 1:7
    plot(age, values[i,:], color = "red", alpha = 0.5, linewidth = 3)
end
for i in 8:14
    plot(age, values[i,:], color = "blue", alpha = 0.5,linewidth = 3)
end

plot(age, surface_values, color = "black", linewidth = 3) 

legend(vcat(num, ["surface"]), loc="center left", bbox_to_anchor=(1, 0.5))
xlabel("Time [yr]")
ylabel("θ anomaly")

using Revise, TMI, GoogleDrive
using PyPlot, PyCall, Test, NCDatasets

filepath = "/home/brynn/Documents/JP/rdata/LGMR_SST_climo.nc"
nc = NCDataset(filepath)
lat = getindex(nc, "lat")
lon = getindex(nc,"lon") #°E
age = getindex(nc,"age") #years BP
sst = nc["sst"][:]

pygui(true)
figure()
contourf(lon[:,1], lat[1,:], sst[:,:,1]')
#issue with tripolar 

#=
Generate similarity matrix 
=#


using TMI
using Revise
using PyPlot, PyCall, CSV, DataFrames, NaNMath

filepath = "/home/brynn/Documents/JP/rdata/coreloc1.csv"
x = CSV.read(filepath,DataFrame,header = 1)
lats = x[:,2]
lons = x[:,3]

#convert to degrees East
lons = 360 .- lons
depths = x[:,4]
nums = x[:,1]

TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

loc = (lons[1],lats[1],depths[1])
origin = surfaceorigin(loc,Alu,γ)

origins = zeros((length(num), length(γ.lon), length(γ.lat)))
for i in 1:length(num)
    xlon = lons[i]
    xlat = lats[i]
    xdepth = depths[i]
    loc = (xlon,xlat,xdepth)
    origins[i, :, :]  .= surfaceorigin(loc, Alu, γ)
end

#compute mse 
function mse(array1, array2)
    n = length(vec(array1))
    ans = NaNMath.sum((log10.(vec(array1)) .- log10.(vec(array2))).^2) ./ n
    return ans
end

mses =  zeros((length(num), length(num)))
for i in 1:length(num)
    for j in 1:length(num)
        m = mse(origins[i,:,:], origins[j,:,:])
        m = m == 0 ? NaN : m 
        mses[i,j] = m 
    end
end

ms = matshow(mses, interpolation="nearest",cmap = "plasma")
xticks(0:(length(num))-1,labels = [string(i) for i in num],  rotation=90)
yticks(0:(length(num))-1,labels = [string(i) for i in num])
grid()
colorbar(ms)

g1_ind = 1:7
g2_ind = 8:length(num)

m1 = log10.( mean(origins[g1_ind,:,:],dims = 1)[1,:,:])

m2 = mean(origins[g2_ind,:,:],dims = 1)[1,:,:]

figure()
subplot(2,1,1)
contourf(γ.lon, γ.lat, m1')
subplot(2,1,2)
contourf(γ.lon,γ.lat,m2')

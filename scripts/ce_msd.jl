#=
Generate similarity matrix
For a list of core locations (in `filepath`)
  (1) generate surface origin density map 
  (2) vectorize 
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

origins = zeros((length(nums), length(γ.lon), length(γ.lat)))
for i in 1:length(nums)
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

mses =  zeros((length(nums), length(nums)))
for i in 1:length(nums)
    for j in 1:length(nums)
        m = mse(origins[i,:,:], origins[j,:,:])
        m = m == 0 ? NaN : m 
        mses[i,j] = m 
    end
end

ms = matshow(mses, interpolation="nearest",cmap = "plasma")
xticks(0:(length(nums))-1,labels = [string(i) for i in nums],  rotation=90)
yticks(0:(length(nums))-1,labels = [string(i) for i in nums])
grid()
colorbar(ms)

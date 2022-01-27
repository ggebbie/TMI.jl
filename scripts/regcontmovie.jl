#=
Understand how a specific sea/region propagates, application of ex2 
  (1) For a region, generate water mass distribution @ steady state 
  (2a) make an animation of how much water is propagated at each depth 
  (2b) 3d Volume map of water mass - still needs work

Note: Makie can't display both these things at once, mp4 will display once 
and is saved within /scripts. Select and run the volume code to generate 
an interactive Makie popup 
=#
using Revise
using TMI, GLMakie 

#load in TMI 
TMIversion = "modern_180x90x33_GH11_GH12"
#TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

#regions of interest for this problem, pick one below 
list = ("NATL","GIN","LAB")
region = list[3]

#compute watermass distribution for region 
g = 100 * watermassdistribution(TMIversion, Alu, region, γ)
g = cat(g[91:end,:,:],g[begin:90, :,:], dims = 1) #center Atlantic (janky!)

#get x, y, z values 
depth = 0
ind = findall(x->x==depth, γ.depth)[1]
lon = -179:2:179
lat = γ.lat
depth = γ.depth

lims = 0:5:105

#make animation code 
#fname = plotsdir()*"/" *region*"_depth.mp4"
fname = region * "_depth.mp4"
fig, ax, cf = GLMakie.contourf(lon, lat, g[:,:,1], colormap = :plasma, levels = lims)
c = GLMakie.contour!(lon, lat, g[:,:,1],levels = lims,color="black")
scatter!([-22],[62], color = :yellow)
GLMakie.Colorbar(fig[1,2], cf)
nframes = length(depth)
framerate = 5
iterator = 1:nframes
ax.title = "Depth = 0"
ax.xlabel = "Latitude [°N]"
ax.ylabel = "Longitude [°E]"
ylims!(ax, 55,80)
xlims!(ax, -60,0)

record(fig, fname, iterator; framerate=framerate) do d
    cf[3] = g[:,:,d]
    c[3] = g[:,:,d]
    ax.title = "Depth = " *string(depth[d])
end

#now create a 3d model - still something weird with this - mirroring? 
g_new = similar(g)
for i in 1:length(lon)
    for j in 1:length(lat)
        for k in 1:length(depth)
            val = g[i,j,k] / 100
            val = val < 0.05 ? NaN : val
            g_new[i,j,k] = val
        end
    end
end

colormap = RGBAf.(to_colormap(:plasma),1)
colormap[1] = RGBAf(0,0,0,0)
volume(lon, lat, depth ./ 100, g_new,algorithm = :absorption,absorption = 4f0, colormap=colormap)

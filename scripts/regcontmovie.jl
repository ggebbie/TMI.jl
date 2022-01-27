#=
depth varying movie of ex2 for cont of LAB vs GIN 
=#
using Revise
using TMI, PyPlot

# pygui(true) #needed for Atom, not sure what it will do in other places

TMIversion = "modern_180x90x33_GH11_GH12"
#TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# region = list[12]
list = ("NATL","GIN","LAB")

region = list[3]
region = 
g = 100 * watermassdistribution(TMIversion, Alu, region, γ)
g = cat(g[91:end,:,:],g[begin:90, :,:], dims = 1)


depth = 0
ind = findall(x->x==depth, γ.depth)[1]

lon = -179:2:179
lat = γ.lat
lims = 0:5:105
depth = γ.depth

#fig, ax, cf = GLMakie.contourf(lon, lat, g[:,:,1], colormap = :plasma, levels = lims)
#c = GLMakie.contour!(lon, lat, g[:,:,1],levels = lims,color="black")
#scatter!([-22],[62], color = :yellow)
#GLMakie.Colorbar(fig[1,2], cf)
#nframes = length(depth)
#framerate = 5
#iterator = 1:nframes
#ax.title = "Depth = 0"
#ax.xlabel = "Latitude [°N]"
#ax.ylabel = "Longitude [°E]"
#ylims!(ax, 55,80)
#xlims!(ax, -60,0)
#record(fig, region*"_depth.mp4", iterator; framerate=framerate) do d
#    cf[3] = g[:,:,d]
#    c[3] = g[:,:,d]
#    ax.title = "Depth = " *string(depth[d])
#end

aind = findall(x->x ∈ (61,63,59), lat)
oind = findall(x->x ∈ (-21,-23,-25), lon)   
gin = 100*watermassdistribution(TMIversion, Alu, "GIN", γ)
gin = cat(gin[91:end,:,:],g[begin:90, :,:], dims = 1)
lab = 100*watermassdistribution(TMIversion, Alu, "LAB", γ)
lab = cat(lab[91:end,:,:],g[begin:90, :,:], dims = 1)

pairs = [(75,78),(76,79),(77,80)]
figure(figsize = (2, 10))
for i in 1:length(pairs)
    subplot(length(pairs),1,i)
    plot(depth, gin[pairs[i][2], pairs[i][1], :])
    plot(depth, lab[pairs[i][2], pairs[i][1], :])
    legend(["GIN", "LAB"])
    if i == 3
        xlabel("Depth [m]")
    end
    
    ylabel("%")
    title(string((lat[pairs[i][1]], lon[pairs[i][2]])))
end

region = "GIN"
g = watermassdistribution(TMIversion, Alu, region, γ)
g = cat(g[91:end,:,:],g[begin:90, :,:], dims = 1)

g_new = similar(g)
for i in 1:length(lon)
    for j in 1:length(lat)
        for k in 1:length(depth)
            val = g[i,j,k]
            val = val < 0.05 ? NaN : val
            g_new[i,j,k] = val
        end
    end
end

colormap = RGBAf.(to_colormap(:plasma),1)
colormap[1] = RGBAf(0,0,0,0)
volume(g_new,algorithm = :absorption,absorption = 4f0, colormap=colormap)

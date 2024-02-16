#=
% Example: Compute the spatial distribution of δ¹³C %
=#
using Revise
using TMI, PyPlot, PyCall, Statistics, DrWatson

mkpath(plotsdir())

TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
δ¹³C = getVar(TMIversion,"δ¹³C")
PO₄ᴿ = getVar(TMIversion,"PO₄")
depth = getVar(TMIversion, "depth")
d_NATL = getRegion(TMIversion,"d_NATL")

# select the surface level (lvl = 1)
lvl = 1

#select a region
d_NATL = getRegion(TMIversion,"d_NATL")

#plot the distribution at the surface
fig = figure()
ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())
cax = ax.contourf(γ.lon,γ.lat,δ¹³C[:, :,lvl]' .* d_NATL' ,
            transform=TMI.cartopy.crs.PlateCarree()
            ) 
ax.coastlines() #show coastlines
gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
gl.top_labels = false
gl.right_labels = false
ax.set_title("δ¹³C d_NATL")
fig.colorbar(cax, orientation = "horizontal")
fig
fig.savefig(plotsdir("δ¹³C_d_NATL.png"))


d_TROPPAC = getRegion(TMIversion,"d_TROPPAC")

fig = figure()
ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())
ax.contourf(γ.lon,γ.lat,δ¹³C[:, :,lvl]' .* d_TROPPAC' ,
            transform=TMI.cartopy.crs.PlateCarree()
            ) 
ax.coastlines() #show coastlines
gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
gl.top_labels = false
gl.right_labels = false
ax.set_title("δ¹³C d_TROPPAC")
fig.colorbar(cax, orientation = "horizontal")
fig
fig.savefig(plotsdir("δ¹³C_d_TROPPAC.png"))


# δ¹³C_AS = δ¹³C_P - δ¹³C_NA

# println("δ¹³C_AS ", δ¹³C_AS)

# x = PO₄ᴿ[:, :,lvl]' .* d_NATL'
# PO₄ᴿ_NA = mean(filter(!isnan,x))
# R = -3 / 2.5

# pred_δ¹³C = ones(length(axes(δ¹³C,3)))
# true_δ¹³C = similar(pred_δ¹³C)
# for level in size(δ¹³C,3)
#     x = PO₄ᴿ[:, :,level]' .* d_SUBANTPAC'
#     PO₄ᴿ_DP = mean(filter(!isnan,x))
#     RΔPO₄ᴿ = R * (PO₄ᴿ_DP - PO₄ᴿ_NA)
#     pred_δ¹³C[level] = δ¹³C_NA + RΔPO₄ᴿ + δ¹³C_AS
#     x = δ¹³C[:, :,level]' .* d_SUBANTPAC'
#     true_δ¹³C[level] = mean(filter(!isnan,x))
# end

# figs, axs = subplots(1,3, figsize = (10, 5))
# axs[1].plot(pred_δ¹³C, -depth)
# axs[1].set_title("pred_δ¹³C")
# axs[2].plot(true_δ¹³C, -depth)
# axs[2].set_title("true_δ¹³C")
# axs[3].plot(pred_δ¹³C .- true_δ¹³C, -depth)
# axs[3].set_title("pred_δ¹³C .- true_δ¹³C")

# println("Compared to: ", mean(filter(!isnan,x)))

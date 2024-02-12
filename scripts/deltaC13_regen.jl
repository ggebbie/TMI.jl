#=
% Example: Compute the spatial distribution of δ¹³C %
=#

using Revise
using TMI, PyPlot, PyCall, Statistics

pygui(true) #needed for Atom, not sure what it will do in other places

# A, Alu, γ, inputfile = config(url,inputdir)
# ΔPO₄ = readtracer(inputfile,"qpo4")
TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
δ¹³C = getVar(TMIversion,"δ¹³C")
PO₄ᴿ = getVar(TMIversion,"PO₄")
depth = getVar(TMIversion, "depth")
d_NATL = getRegion(TMIversion,"d_NATL")
# plot a section at 330 east longitude (i.e., 30 west)

# view the surface
# cntrs = 1:0.25:6

# PyPlot turned off for CI.
# last index is the depth level? so the is complete and we can compute as 
# planned
lvl = 1

d_NATL = getRegion(TMIversion,"d_NATL")
fig = figure()
ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())
ax.contourf(γ.lon,γ.lat,δ¹³C[:, :,lvl]' .* d_NATL' ,
            transform=TMI.cartopy.crs.PlateCarree()
            ) 
ax.coastlines() #show coastlines
gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
gl.top_labels = false
gl.right_labels = false
ax.set_title("δ¹³C d_NATL")

x = δ¹³C[:, :,lvl]' .* d_NATL'
δ¹³C_NA = mean(filter(!isnan,x)) 
println("δ¹³C_NA: ", δ¹³C_NA)
#d_TROPPAC
# d_SUBANT= getRegion(TMIversion,"d_SUBANT")
d_SUBANTPAC = getRegion(TMIversion,"d_SUBANTPAC")
fig = figure()
ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())
ax.contourf(γ.lon,γ.lat,δ¹³C[:, :,lvl]' .* d_SUBANTPAC' ,
            transform=TMI.cartopy.crs.PlateCarree()
            ) 
ax.coastlines() #show coastlines
gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
gl.top_labels = false
gl.right_labels = false
ax.set_title("δ¹³C d_TROPPAC")
x = δ¹³C[:, :,lvl]' .* d_SUBANTPAC'
δ¹³C_P = mean(filter(!isnan,x))
println("δ¹³C_P: ", δ¹³C_P)

δ¹³C_AS = δ¹³C_P - δ¹³C_NA

println("δ¹³C_AS ", δ¹³C_AS)

x = PO₄ᴿ[:, :,lvl]' .* d_NATL'
PO₄ᴿ_NA = mean(filter(!isnan,x))
R = -3 / 2.5

pred_δ¹³C = ones(length(axes(δ¹³C,3)))
true_δ¹³C = similar(pred_δ¹³C)
for deep_lvl in axes(δ¹³C,3)
    x = PO₄ᴿ[:, :,deep_lvl]' .* d_SUBANTPAC'
    PO₄ᴿ_DP = mean(filter(!isnan,x))
    RΔPO₄ᴿ = R * (PO₄ᴿ_DP - PO₄ᴿ_NA)
    pred_δ¹³C[deep_lvl] = δ¹³C_NA + RΔPO₄ᴿ + δ¹³C_AS
    x = δ¹³C[:, :,deep_lvl]' .* d_SUBANTPAC'
    true_δ¹³C[deep_lvl] = mean(filter(!isnan,x))
end

figs, axs = subplots(1,3, figsize = (10, 5))
axs[1].plot(pred_δ¹³C, -depth)
axs[1].set_title("pred_δ¹³C")
axs[2].plot(true_δ¹³C, -depth)
axs[2].set_title("true_δ¹³C")
axs[3].plot(pred_δ¹³C .- true_δ¹³C, -depth)
axs[3].set_title("pred_δ¹³C .- true_δ¹³C")

println("Compared to: ", mean(filter(!isnan,x)))
# fig = figure()
# ax = fig.add_subplot(projection = TMI.cartopy.crs.PlateCarree())
# ax.contourf(γ.lon,γ.lat,PO₄ᴿ[:, :, lvl]') # units: effective thickness in log10(meters)
# ax.coastlines() #show coastlines
# gl = ax.gridlines(draw_labels=true, dms=true, x_inline=false, y_inline=false)
# gl.top_labels = false
# gl.right_labels = false
# ax.set_title("PO₄ᴿ")


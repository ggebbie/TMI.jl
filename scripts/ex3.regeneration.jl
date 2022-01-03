#=
% Example: Determine the total amount of remineralized phosphate   %
=#

using Revise
using TMI, PyPlot, PyCall

pygui(true) #needed for Atom, not sure what it will do in other places

# A, Alu, γ, inputfile = config(url,inputdir)
# ΔPO₄ = readtracer(inputfile,"qpo4")
TMIversion = "modern_90x45x33_GH10_GH12"
#TMIversion = "modern_180x90x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
PO₄ᴿ = regeneratedphosphate(TMIversion,Alu,γ)


# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330; # only works if exact
Psection = section(PO₄ᴿ,lon_section,γ)
lims = 0:0.1:2.0

# make a plot of regenerated phosphate
dyeplot(γ.lat,-γ.depth[33:-1:1],Psection[:,33:-1:1]', lims)

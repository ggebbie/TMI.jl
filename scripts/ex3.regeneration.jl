#=
% Example: Determine the total amount of remineralized phosphate   %
=#

import Pkg; Pkg.activate(".")

using Revise
using TMI
using Test
using GeoPythonPlot
#using BenchmarkTools

# A, Alu, γ, inputfile = config(url,inputdir)
# ΔPO₄ = readtracer(inputfile,"qpo4")
TMIversion = "modern_90x45x33_GH10_GH12"
#TMIversion = "modern_180x90x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

PO₄ᴿ = regeneratedphosphate(TMIversion,Alu,γ)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330 # only works if exact
lims = 0:0.1:2.0
tlabel = "Regenerated Phosphate [μmol/kg]"
sectionplot(PO₄ᴿ,lon_section,lims,titlelabel = tlabel)

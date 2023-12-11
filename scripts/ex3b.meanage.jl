#=
% Diagnose the mean or ideal age.
=#

import Pkg; Pkg.activate(".")

using Revise
using TMI
using Test
using GeoPythonPlot

#TMIversion = "modern_90x45x33_GH10_GH12"
TMIversion = "modern_180x90x33_GH11_GH12"

A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

a = meanage(TMIversion,Alu,γ)

# plot a section at 330 east longitude (i.e., 30 west)
#lon_section = 330 # only works if exact
#lon_section = -150
lon_section = -151
lims = 0:100:2000
tlabel = string(lon_section)*"°E, Mean Age [yr]"
sectionplot(a,lon_section,lims,titlelabel = tlabel)

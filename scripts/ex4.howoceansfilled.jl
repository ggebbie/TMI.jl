#=
% Example: Find the ocean volume that has originated from each    %
%            surface box.                                           %
%                                                                   %
% This is equivalent to solving a sensitivity problem:              %
% The total volume is V = v^T c , where v is the volume of each box %
% and c is the fraction of volume from a given source which         %
% satisfies the equation A c = d.                                   %
% Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
% originating from each source.
%
% See Section 3 and Supplementary Section 4, Gebbie & Huybers 2011. 
=#

import Pkg; Pkg.activate("./Documents/GitHub/TMI.jl/scripts")

using Revise
using TMI
using Test
# using GeoPythonPlot
using PythonCall
using CondaPkg
using PythonPlot
const cartopy = pyimport("cartopy")
ccrs = cartopy.crs

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);


volume = volumefilled(TMIversion,Alu,γ)

# view the surface
cntrs = 1:0.25:6
tlabel="Volume filled by surface patch [log₁₀(m)]"
GeoPythonPlot.pygui(true) # to help plots appear on screen using Python GUI
planviewplot(volume,cntrs,titlelabel=tlabel) 

fig, axs = subplots(subplot_kw=Dict("projection"=> ccrs.PlateCarree()), figsize=(11,8.5))


axs.coastlines()
fig














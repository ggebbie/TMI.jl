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

import Pkg; Pkg.activate("./scripts")
Pkg.instantiate()
# using Revise # this package is not included in scripts toml files not added in 
using TMI
using Test
using GeoPythonPlot

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

volume = volumefilled(TMIversion,Alu,γ)

lon = volume.axes[1]
lat = volume.axes[2]
z = 0 #we are looking at surface sources
# view the surface
cntrs = 1:0.25:6
tlabel="Volume filled by surface patch [log₁₀(m)]"

GeoPythonPlot.pygui(true) # to help plots appear on screen using Python GUI
planviewplot(volume.tracer, lon, lat, z,cntrs,titlelabel=tlabel) 

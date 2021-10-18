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
using Revise
using TMI, PyPlot, PyCall

TMIversion = "TMI_2010_2012_4x4x33"
volume = volumefilled(TMIversion)

# view the surface
cntrs = 1:0.25:6

# PyPlot turned off for CI.
contourf(γ.lon,γ.lat,log10.(volume'),cntrs) # units: effective thickness in log10(meters)

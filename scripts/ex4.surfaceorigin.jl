#= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Example: Find the surface origin of water for some interior box %
 %                                                                   %
 % This is equivalent to solving a sensitivity problem:              %
 % The total volume is V = v^T c , where v is the volume of a
 % given interior box,
 % and c is the fraction of volume from a given source which         %
 % satisfies the equation A c = δ.                                   %
 % Next, dV/d(δ) = A^(-T) v, and dV/d(δ) is exactly the volume       %
 % originating from each source.      
 % Very similar mathematically to determining how the ocean is filled.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
=#
using Revise
using TMI, PyPlot, PyCall

# % choose an interior location X (Xlon[lon], Xlat [lat], Xdepth [m depth]).
#% -7.38, 115.26E

xlon = 125.26; # deg E.
xlat = -6.38;  # deg N.
xdepth = 3000;  # meters.
loc = (xlon,xlat,xdepth)

origin, γ = surfaceorigin(TMIversion,loc)

# units: effective thickness in log10(meters)
contourf(γ.lon,γ.lat,log10.(origin'))

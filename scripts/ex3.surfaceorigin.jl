#= %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Example 3: Find the surface origin of water for some interior box %
 %                                                                   %
 % This is equivalent to solving a sensitivity problem:              %
 % The total volume is V = v^T c , where v is the volume of a
 % given interior box,
 % and c is the fraction of volume from a given source which         %
 % satisfies the equation A c = d.                                   %
 % Next, dV/d(d) = A^(-T) v, and dV/d(d) is exactly the volume       %
 % originating from each source.      
 % Very similar mathematically to example 2.                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
=#
using Revise
using TMI, Interpolations, PyPlot, PyCall

# % choose an interior location X (Xlon[lon], Xlat [lat], Xdepth [m depth]).
#% -7.38, 115.26E

xlon = 125.26; # deg E.
xlat = -6.38;  # deg N.
xdepth = 3000;  # meters.
loc = (xlon,xlat,xdepth)

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../data"

A, Alu, c, γ = config(url,inputdir)

# Find nearest neighbor on grid
# set δ = 1 at grid cell of interest
δ = nearestneighbormask(loc,γ)

# Include option: find grid coordinate by linear interpolation/extrapolation
# try Interpolations.jl, need interpolation factors that add up to one
# a false start to fix this issue 
#v = cellVolume(γ)
#itp = interpolate(v)

#
dVdδ = tracerinit(γ.wet); # pre-allocate c
dVdδ[γ.wet] = Alu'\δ[γ.wet]

# surfaceorigin only exists at sea surface
surfaceorigin = view(dVdδ,:,:,1)

# view the surface
#cntrs = 1:0.25:6

# units: effective thickness in log10(meters)
figure()
contourf(γ.lon,γ.lat,log10.(surfaceorigin'))

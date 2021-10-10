#=
% Example 2: Find the ocean volume that has originated from each    %
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
using TMI, PyPlot, PyCall, BenchmarkTools,
 GibbsSeaWater, Distances 

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../data"

Azyx, Axyz, Alu, c, γ = config(url,inputdir)

v = cellVolume(γ)
area = cellArea(γ)
 
# effectively take inverse of transpose A matrix.
dVdd = tracerFieldInit(γ.wet); # pre-allocate c
dVdd[γ.wet] = Alu'\v[γ.wet]

# scale the sensitivity value by surface area so that converging meridians are taken into account.
I = γ.I
volumefill = Matrix{Float64}(undef,length(γ.lon),length(γ.lat))
fill!(volumefill,0.0)

# this step could use a function with γ.I argument
[volumefill[I[ii][1],I[ii][2]] = dVdd[I[ii]] / area[I[ii][1],I[ii][2]] for ii ∈ eachindex(I) if I[ii][3] == 1]

# view the surface
cntrs = 1:0.25:6

# PyPlot turned off for CI.
contourf(γ.lon,γ.lat,log10.(volumefill'),cntrs) # units: effective thickness in log10(meters)

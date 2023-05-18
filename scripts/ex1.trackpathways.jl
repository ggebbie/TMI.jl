#=
 Example 1: Track the pathways of a user-defined water mass.
 Steps: (a) define the water mass 1). by a rectangular surface patch
            dyed with passive tracer concentration of 1,
            or 2. load a pre-defined surface patch in d_all.mat.
        (b) propagate the dye with the matrix A, with the result
            being the fraction of water originating from the
            surface region.
 See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
=#

import Pkg; Pkg.activate(".")

using Revise
using TMI
using Test
using GGplot

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

#- define the surface patch by the bounding latitude and longitude.
latbox = [50,60]; # 50 N -> 60 N, for example.

# mutable due to wraparound: don't use an immutable tuple
lonbox = [-50,0]; # 50 W -> prime meridian

# do numerical analysis
c = trackpathways(Alu,latbox,lonbox,γ);

# do plotting (could be a function)
plotextent(latbox, lonbox)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330; # only works if exact
lims = 0:5:100
tlabel = "water-mass fraction [%]"
sectionplot(100c,lon_section,lims,titlelabel = tlabel)

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
using Revise
using TMI, BenchmarkTools, PyPlot, PyCall

pygui(true) #needed for Atom, not sure what it will do in other places

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../data"

#c = readTracer(url,"θ")

A, Alu, c, ΔPO₄, γ = config(url,inputdir)

#- define the surface patch by the bounding latitude and longitude.
latbox = [50,60]; # 50 N -> 60 N, for example.

# mutable due to wraparound: don't use an immutable tuple
lonbox = [-50,0]; # 50 W -> prime meridian

d = surfacepatch(lonbox,latbox,γ)

# do matrix inversion to get quantity of dyed water throughout ocean:
c = tracerinit(γ.wet); # pre-allocate c

# make methods that make the "wet" index unnecessary
c[γ.wet] = Alu\d[γ.wet] # presumably equivalent but faster than `c = A\d`

#plot bbox
plotextent(latbox, lonbox)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330;
csection = section(c,lon_section,γ)
lims = 0:5:100

# make a plot of dye in the ocean
dyeplot(γ.lat,-γ.depth[33:-1:1],100 * csection[:,33:-1:1]', lims)

using Revise, TMI
using Plots
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

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../input"

A, Alu, γ = configTMI(url,inputdir)

#- define the surface patch by the bounding latitude and longitude.
latbox = (50,60); # 50 N -> 60 N, for example.

# mutable due to wraparound: don't use an immutable tuple
lonbox = [-50,0]; # 50 W -> prime meridian

# ternary operator to handle longitudinal wraparound
lonbox[1] ≤ 0 ? lonbox[1] += 360 : nothing
lonbox[2] ≤ 0 ? lonbox[2] += 360 : nothing

d = surfacepatch(lonbox,latbox,γ)

# do matrix inversion to get quantity of dyed water throughout ocean:
c = Alu\d

# after doing calculations with vectors, translate to a 3D geometric field
# must be easier way to do this with cartesianindex type
cfld = vec2fld(c,γ.I)

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330;
isec = findall(==(lon_section),γ.lon)
csection = dropdims(cfld[isec,:,:],dims=1)

# make a plot of dye in the ocean
Plots.contourf(γ.lat,-γ.depth[33:-1:1],csection[:,33:-1:1]') # a sample plot at 22 W.

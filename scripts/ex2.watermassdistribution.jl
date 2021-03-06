#=
 Example 2: Find the global distribution of a "water mass" defined by an oceanographically-relevant surface region.
 Steps: (a) define the water mass 1). by a pre-defined surface
            dyed with passive tracer concentration of 1,
        (b) propagate the dye with the matrix A, with the result
            being the fraction of water originating from the
            surface region.
 See Section 2b of Gebbie & Huybers 2010, esp. eqs. (15)-(17).
=#
using Revise
using TMI, BenchmarkTools, PyPlot, PyCall

#TMIversion = "modern_180x90x33_GH10_GH12"
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# choose water mass (i.e., surface patch) of interest
# Enhancement: provide list of choices
list = ("GLOBAL","ANT","SUBANT",
            "NATL","NPAC","TROP","ARC",
            "MED","ROSS","WED","LAB","GIN",
            "ADEL","SUBANTATL","SUBANTPAC","SUBANTIND",
            "TROPATL","TROPPAC","TROPIND")

region = list[2]

# do numerical analysis
g = watermassdistribution(TMIversion,Alu,region,γ);

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330 # only works if exact
lims = 0:5:100
tlabel = region * " water-mass fraction [%]"
sectionplot(100g,lon_section,lims,titlelabel = tlabel)

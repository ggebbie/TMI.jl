#=
% Example: Determine the total amount of remineralized phosphate   %
=#

using Revise
using TMI, PyPlot, PyCall

pygui(true) #needed for Atom, not sure what it will do in other places

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../data"

#c = readTracer(url,"θ")

A, Alu, c, ΔPO₄, γ = config(url,inputdir)

# PO₄ᴿ = cumulative regenerated phosphate
PO₄ᴿ = tracerinit(γ.wet); # pre-allocate 
PO₄ᴿ[γ.wet] = -(Alu\ΔPO₄[γ.wet])

# plot a section at 330 east longitude (i.e., 30 west)
lon_section = 330;
Psection = section(PO₄ᴿ,lon_section,γ)
lims = 0:0.1:2.0

# make a plot of regenerated phosphate
dyeplot(γ.lat,-γ.depth[33:-1:1],Psection[:,33:-1:1]', lims)

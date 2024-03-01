import Pkg; Pkg.activate(".")
using Revise
using TMI
fname = TMI.pkgdatadir("simple_MITgcm.nc")
println(fname)
γ =  Grid(fname, "maskC", "XC", "YC", "Z")
y = (
    θ = readfield(fname,"THETA",γ),
    S = readfield(fname,"SALT",γ))
#TMI.meanage 
@time m̃ = massfractions(y)
massfractions(y)

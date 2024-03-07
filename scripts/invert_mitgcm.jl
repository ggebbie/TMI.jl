import Pkg; Pkg.activate(".")
using Revise
using TMI
fname = TMI.pkgdatadir("simple_MITgcm.nc")
println(fname)
γ =  Grid(fname, "maskC", "XC", "YC", "Z")

y = (θ = readfield(fname,"THETA",γ),
    S = readfield(fname,"SALT",γ))
    
w = (θ =  0.01,
    S = 0.001)

@time m̃ = massfractions(y, w)
#massfractions(y)

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
Ã = watermassmatrix(m̃, γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
bθ = getsurfaceboundary(y.θ)

## reconstruct temperature
using LinearAlgebra: lu
Ãlu = lu(Ã)
θ̃ = steadyinversion(Ãlu,bθ,γ)

# compare to c.θ
Base.maximum(y.θ - θ̃)
Base.minimum(y.θ - θ̃)


#massfractions(y)

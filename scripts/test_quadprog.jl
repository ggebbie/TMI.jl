using TMI: massfractions_isotropic
import Pkg; Pkg.activate(".")
using Revise
using TMI
fname = TMI.pkgdatadir("simple_MITgcm.nc")
println(fname)
γ =  Grid(fname, "maskC", "XC", "YC", "Z")
y = (θ = readfield(fname,"THETA",γ),
    S = readfield(fname,"SALT",γ))

# m̃ = massfractions(y)

# pull this function apart
c = y 
γ = first(c).γ
Rfull = CartesianIndices(γ.wet)
Ifirst, Ilast = first(Rfull), last(Rfull)

    # loop over interior points (dirichlet bc at surface)
    R = copy(γ.R)
    Iint = cartesianindex(γ.interior)

    m0 = TMI.massfractions_isotropic(γ)
    @time n   = TMI.neighbors(m0,γ)
    @time n2 = TMI.neighbors(γ)
    nrow   = length(c) + 1 # add mass conservation

    # allocate maximum needed
    nmax = maximum(neighbors)
    l = zeros(nmax)
    u = ones(nmax)
    b = vcat(zeros(nrow-1),1.0)
    mlocal = zeros(nmax)
    nlocal = zeros(nrow)
    ϵ = 1e-8 # for checking tolerances



massfractions


import Pkg; Pkg.activate(".")

using TMI
using Revise
using LinearAlgebra, Statistics, SparseArrays
# using GeoPythonPlot


versions = ["LGM_90x45x33_G14", "LGM_90x45x33_OG18", "LGM_90x45x33_GPLS2", "LGM_90x45x33_G14A"]

for TMIversion in versions
    println(TMIversion)
    TMIversion_compressed =  TMIversion*"_compressed.nc"
    TMIfile_compressed = TMI.pkgdatadir("TMI_" * TMIversion_compressed)
    γ = Grid(TMIfile_compressed)
    θ = readfield(TMIfile_compressed,"θ",γ)
    θ_ref = readfield(TMIfile_compressed,"θ_reference",γ)
    println(mean(θ))
    println(mean(θ_ref))
end

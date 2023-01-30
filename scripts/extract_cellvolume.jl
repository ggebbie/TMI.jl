#=
% Extract TMI grid volume for ocean points in 2x2 and 4x4 resolutions.
=#

import Pkg; Pkg.activate(".")

using TMI
using MAT

versions = ("modern_90x45x33_GH10_GH12","modern_180x90x33_GH11_GH12")

volume = Dict{String,Array{Float64,3}}()
for vers in versions
    A, Alu, γ, TMIfile, L, B = config_from_nc(vers)
    v = cellvolume(γ)

    # save to a Dictionary
    volume[vers] = v
end

fname = "TMI_volumes.mat"
matwrite("../data/"*fname, volume)  

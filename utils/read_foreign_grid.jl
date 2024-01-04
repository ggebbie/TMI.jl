#=
% Create a TMI.Grid struct from foreign model output
=#

import Pkg; Pkg.activate("../scripts") # a more full-featured environment

using Revise
using TMI
using NCDatasets


gamma =  Grid(TMI.pkgdatadir("simple_MITgcm.nc"), "THETA", "XC", "YC", "Z", "maskC") 

# scratch work: shows how the Grid constructor works
mitgcm_grid = Dataset(TMI.pkgdatadir("mitgcm_grid.nc"))
simple_mitgcm = Dataset(TMI.pkgdatadir("simple_MITgcm.nc"))

# use potential temperature as a template
# should double-check consistency with other tracers (but is not done here)
# problem: Union Missing type
c = convert(Array{Float32,3},simple_mitgcm["THETA"][:,:,:,begin]) # take one snapshot
lon = convert(Vector{Float32},simple_mitgcm["XC"]) 
lat = convert(Vector{Float32},simple_mitgcm["YC"])
depth = - convert(Vector{Float32},simple_mitgcm["Z"]) # flip sign for actual "depth"
wet = Bool.(simple_mitgcm["maskC"])::BitArray # very slow! (couple of secs), use `convert` instead?

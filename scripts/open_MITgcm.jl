#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Under development: the goal is to invert model output 
%                  and get the transport matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#

import Pkg; Pkg.activate(".")
using Revise
using TMI
using Test
using GeoPythonPlot
using LinearAlgebra
using SparseArrays
using Statistics
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "simple_mitgcm"

  
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

TMIfile = pkgdatadir(""*TMIversion*".nc")
println(TMIfile)
γ = Grid(TMIfile, "THETA", "XC", "YC", "Z", "maskC")
ds = Dataset(TMIfile)
#θtrue = convert(Vector{Float32},ds["THETA"]) 
θtrue = readfield(TMIfile,"THETA",γ)

#println(γ)

# plot the difference
level = 15 # your choice 1-33
depth = γ.depth[level]

cntrs = 0:0.5:15
label = "True θ"
planviewplot(θtrue, depth, cntrs, titlelabel=label)
readline()

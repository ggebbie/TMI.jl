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
using JLD
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "simple_mitgcm"

  
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

TMIfile = pkgdatadir(""*TMIversion*".nc")
println(TMIfile)
γ = Grid(TMIfile, "hFacC", "XC", "YC", "Z")
θtrue = readfield(TMIfile,"THETA",γ)

i_flatten = lonindex(γ.I)
j_flatten = latindex(γ.I)
k_flatten = depthindex(γ.I)

#println(γ.I)
ind_1 = zeros(0)
ind_2 = zeros(0)

# initialize A
for ii in 1:length(i_flatten)
      loc = zeros(0)
      loc1 = findall((lonindex(γ.I) .== i_flatten[ii]) .& (latindex(γ.I) .== j_flatten[ii]) .& (depthindex(γ.I) .== k_flatten[ii]))
      append!(loc,loc1)
      if i_flatten[ii]>1
          loc1 = findall((lonindex(γ.I) .== i_flatten[ii]-1) .& (latindex(γ.I) .== j_flatten[ii]) .& (depthindex(γ.I) .== k_flatten[ii]))
          append!(loc,loc1)
      end
      if i_flatten[ii]<length(γ.lon)
          loc1 = findall((lonindex(γ.I) .== i_flatten[ii]+1) .& (latindex(γ.I) .== j_flatten[ii]) .& (depthindex(γ.I) .== k_flatten[ii]))
          append!(loc,loc1)
      end
      if i_flatten[ii]==1
          loc1 = findall((lonindex(γ.I) .== length(γ.lon)) .& (latindex(γ.I) .== j_flatten[ii]) .& (depthindex(γ.I) .== k_flatten[ii]))
          append!(loc,loc1)
      end
      if j_flatten[ii]>1
          loc1 = findall((lonindex(γ.I) .== i_flatten[ii]) .& (latindex(γ.I) .== j_flatten[ii]-1) .& (depthindex(γ.I) .== k_flatten[ii]))
          append!(loc,loc1)
      end
      if j_flatten[ii]<length(γ.lat)
          loc1 = findall((lonindex(γ.I) .== i_flatten[ii]) .& (latindex(γ.I) .== j_flatten[ii]+1) .& (depthindex(γ.I) .== k_flatten[ii]))
          append!(loc,loc1)
      end
      if k_flatten[ii]>1
          loc1 = findall((lonindex(γ.I) .== i_flatten[ii]) .& (latindex(γ.I) .== j_flatten[ii]) .& (depthindex(γ.I) .== k_flatten[ii]-1))
          append!(loc,loc1)
      end
      if k_flatten[ii]<length(γ.depth)
          loc1 = findall((lonindex(γ.I) .== i_flatten[ii]) .& (latindex(γ.I) .== j_flatten[ii]) .& (depthindex(γ.I) .== k_flatten[ii]+1))
          append!(loc,loc1)
      end      
      append!(ind_1,ii*ones(length(loc)))
      append!(ind_2,loc)
end

println(size(ind_1),size(ind_2))
save("ind_1.jld","ind_1", ind_1)
save("ind_2.jld","ind_2", ind_2)

# plot the difference
level = 15 # your choice 1-33
depth = γ.depth[level]
lon = γ.lon
lat = γ.lat

cntrs = 0:0.5:15
label = "True θ"
planviewplot(θtrue, depth, cntrs, titlelabel=label,cenlon=40)
readline()

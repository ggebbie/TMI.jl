#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Under development: the goal is to invert model output 
%                  and get the transport matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#

import Pkg; Pkg.activate(".")
using Revise
using TMI
using Test
using GGplot
using LinearAlgebra
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# first guess of change to surface boundary conditions
# ocean values are 0
#u = zerosurfaceboundary(γ)
u = (;surface = zerosurfaceboundary(γ))
uvec = vec(u)
  
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")
θtrue = readfield(TMIfile,"θ",γ) 

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
#b = getsurfaceboundary(y)
b = (;surface = getsurfaceboundary(θtrue))

#We need an error covariance matrix
W⁻ = (1/sum(γ.wet)) .* Diagonal(1 ./(10^(-2) .* ones(sum(γ.wet))).^2)

# get sample J value
F,G = costfunction_gridded_model(uvec,Alu,b,u,θtrue,W⁻,γ)
fg!(F,G,x) = costfunction_gridded_model!(F,G,x,Alu,b,u,θtrue,W⁻,γ)

#fg(x) = costfunction_gridded_model(x,Alu,b,u,θtrue,W⁻,γ)
#f(x) = fg(x)[1]
#J₀,gJ₀ = fg(uvec)

#### gradient check ###################
# check with forward differences
ϵ = 1e-3
ii = rand(1:sum(γ.wet[:,:,1]))
println("Location for test =",ii)
δu = copy(uvec); δu[ii] += ϵ
#∇f_finite = (f(δu) - f(uvec))
#println(∇f_finite)

#fg!(J₀,gJ₀,(uvec+δu)./2) # J̃₀ is not overwritten
#∇f = gJ₀[ii]
#println(∇f)

# error less than 10 percent?
#println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
#### end gradient check #################

# filter the data with an Optim.jl method
iterations = 25
out = steadyclimatology(uvec,fg!,iterations)

# reconstruct by hand to double-check.
ũ = unvec(u,out.minimizer)

# apply to the boundary conditions
b̃ = adjustboundarycondition(b,ũ)

# reconstruct tracer map
c₀ = steadyinversion(Alu,b,γ)
c̃  = steadyinversion(Alu,b̃,γ)

Δc̃ = c̃ - θtrue
Δc₀ = c₀ - θtrue

# plot the difference
level = 15 # your choice 1-33
depth = γ.depth[level]
cntrs = -1:0.05:1
label = "Optimized misfit: Δc̃"
# Help: needs work with continents and labels
planviewplot(Δc̃, depth, cntrs, titlelabel=label) 
readline()

label = "First guess misfit: Δc₀"
planviewplot(Δc₀, depth, cntrs, titlelabel=label) 
readline()

cntrs = -1:0.5:15
label = "True θ"
planviewplot(θtrue, depth, cntrs, titlelabel=label)
readline()

cntrs = -1:0.5:15
label = "Optimized θ"
planviewplot(c̃, depth, cntrs, titlelabel=label)

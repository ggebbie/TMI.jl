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
using SparseArrays
using Statistics
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)


  
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

TMIfile = pkgdatadir("TMI_"*TMIversion*".nc")
θtrue = readfield(TMIfile,"θ",γ)

ctrue = vec(θtrue)

q = -A * ctrue

# The first guess for the tracer concentration should be close to the actual tracer concentration
# take first guess as θtrue+0.01
cvec=vec(θtrue).+ 0.01
θguess = unvec(θtrue,cvec)

#first guess tracer control vector is near zero, and we want this to remain relatively small
u =  Field(-0.01.*ones(size(γ.wet)),γ,θtrue.name,θtrue.longname,θtrue.units)
uvec = vec(u)



#We need an error covariance matrix
W⁻ = Diagonal(1 ./( ones(sum(γ.wet))).^2)#(1/sum(γ.wet))

#I want to allow a bunch of error in the surface part of the tracer conservation please
Qerror = ones(size(γ.wet))
Qerror[:,:,1].=0
Qfield = Field(Qerror,γ,θtrue.name,θtrue.longname,θtrue.units)
Qvec = vec(Qfield)


Q⁻ = Diagonal(1 ./( ones(sum(γ.wet))).^2)

non_zero_indices1, non_zero_indices2, non_zero_values = findnz(A)

non_zero_indices = hcat(non_zero_indices1, non_zero_indices2)


convec = [uvec; non_zero_values-non_zero_values]
A0=A
ulength=length(uvec)

# get sample J value
F = costfunction_gridded_model(convec,non_zero_indices,u,A0,ctrue,cvec,q,W⁻,Q⁻,γ)
fg!(F,G,x) = costfunction_gridded_model!(F,G,x,non_zero_indices,u,A0,ctrue,cvec,q,W⁻,Q⁻,γ)


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


print(length(convec))
# filter the data with an Optim.jl method
iterations = 5
out = steadyclimatology(convec,fg!,iterations)

# reconstruct by hand to double-check.
ũ = unvec((W⁻ * u),out.minimizer)

# apply to the boundary conditions
#b̃ = adjustboundarycondition(b,ũ)

# reconstruct tracer map
c₀ = θguess
c̃  = θguess+ũ

Δc̃ = c̃ - θtrue
Δc₀ = θguess - θtrue

Anew = sparse(non_zero_indices[:, 1], non_zero_indices[:, 2], convec[ulength+1:end])



# plot the difference
level = 15 # your choice 1-33
depth = γ.depth[level]
cntrs = -15:0.5:15
label = "Optimized misfit: Δc̃"
## Help: needs work with continents and labels
planviewplot(Δc̃, depth, cntrs, titlelabel=label) 
readline()

cntrs = -15:0.5:15
label = "First guess misfit: Δc₀"
planviewplot(Δc₀, depth, cntrs, titlelabel=label) 
readline()

cntrs = 0:0.5:15
label = "True θ"
planviewplot(θtrue, depth, cntrs, titlelabel=label)
readline()

cntrs = 0:0.5:15
label = "Optimized θ"
planviewplot(c̃, depth, cntrs, titlelabel=label)



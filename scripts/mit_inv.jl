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

ctrue = Float64.(vec(θtrue))

#println(typeof(ctrue))


surfind = surfaceindex(γ.I)

# The first guess for the tracer concentration should be close to the actual tracer concentration
# take first guess as θtrue
cvec=(vec(θtrue))
θguess = unvec(θtrue,cvec)


#first guess tracer control vector is near zero, and we want this to remain relatively small
u =  Field(-0.01.*ones(size(γ.wet)),γ,θtrue.name,θtrue.longname,θtrue.units)
uvec = vec(u)
q = zeros(size(uvec))
q[surfind] = cvec[surfind]

#We need an error covariance matrix
W⁻ = Diagonal(1 ./( ones(sum(γ.wet))).^2)#(1/sum(γ.wet))


Qdiag = 1 ./ ones(sum(γ.wet))
#Qdiag[surfind] .= 10^(-7)
Q⁻ = Diagonal(Qdiag)

ind_1 = load("ind_1.jld")["ind_1"]
ind_2 = load("ind_2.jld")["ind_2"]

A0 = sparse(ind_1, ind_2, ones(length(ind_1))*0.5)

non_zero_indices1, non_zero_indices2, non_zero_values = findnz(A0)

non_zero_indices = hcat(non_zero_indices1, non_zero_indices2)


convec = [uvec; non_zero_values]
ulength=length(uvec)

# get sample J value
F = costfunction_gridded_model(convec,non_zero_indices,u,A0,ctrue,cvec,q,W⁻,Q⁻,γ)
fg!(F,G,x) = costfunction_gridded_model!(F,G,x,non_zero_indices,u,A0,ctrue,cvec,q,W⁻,Q⁻,γ)
fg(x) = costfunction_gridded_model(x,non_zero_indices,u,A0,ctrue,cvec,q,W⁻,Q⁻,γ)
f(x) = fg(x)[1]
J₀,gJ₀ = fg(convec)

#### gradient check ###################
# check with forward differences
#ϵ = 1e-3
#println(size(length(convec)))
#ii = rand(1:length(convec))
#println("Location for test =",ii)
#δu = copy(convec); δu[ii] += ϵ
#∇f_finite = (f(δu) - f(convec))/ϵ
#println(∇f_finite)

#fg!(J₀,gJ₀,(convec+δu)./2) # J̃₀ is not overwritten
#∇f = gJ₀[ii]
#println(∇f)

# error less than 10 percent?
#println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
#### end gradient check #################

#print(length(convec))
# filter the data with an Optim.jl method
iterations = 2
out = steadyclimatology_optim(convec,fg!,iterations)

# reconstruct by hand to double-check.
ũ = unvec(θguess,Float32.(out.minimizer[begin:ulength]))

#println(typeof(ũ))

# reconstruct tracer map
c₀ = θguess
c̃  = θguess+ũ

Δc̃ = c̃ - θtrue
Δc₀ = θguess - θtrue

Anew = A0 + sparse(non_zero_indices[:, 1], non_zero_indices[:, 2], out.minimizer[ulength+1:end])
#onesvec = ones(size(q))

#oldf = sum((non_zero_values).^2)
#newf = sum((out.minimizer[ulength+1:end]).^2)
#tracer_cons1 = sum((Q⁻*A0*cvec - q).^2)
#tracer_cons2 = sum((Q⁻*Anew*(cvec+out.minimizer[begin:ulength]) - q).^2)
#mass_cons1 = sum((A0*onesvec).^2)
#mass_cons2 = sum((Anew*onesvec).^2)

#println("old tracer cons:$tracer_cons1")
#println("new tracer cons:$tracer_cons2")
#println("old mass cons:$mass_cons1")
#println("new mass cons:$mass_cons2")

b = getsurfaceboundary(θtrue)

#A0


#Anew
#Alu2 = lu(A0)


final_θ = steadyinversion(Anew,b,γ)

# plot the difference
#level = 15 # your choice 1-33
#depth = γ.depth[level]

#cntrs = 0:0.5:15
#label = "True θ"
#planviewplot(θtrue, depth, cntrs, titlelabel=label)
#readline()

#cntrs = 0:0.5:15
#label = "Final θ"
#planviewplot(final_θ, depth, cntrs, titlelabel=label)
#readline()

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

BLAS.set_num_threads(8)


TMIversion = "simple_MITgcm"

  
pkgdir() = dirname(dirname(pathof(TMI)))
pkgdir(args...) = joinpath(pkgdir(), args...)

pkgdatadir() = joinpath(pkgdir(),"data")
pkgdatadir(args...) = joinpath(pkgdatadir(), args...)

TMIfile = pkgdatadir(""*TMIversion*".nc")
println(TMIfile)
γ = Grid(TMIfile, "hFacC", "XC", "YC", "Z")
θtrue = readfield(TMIfile,"THETA",γ)

ctrue = Float64.(vec(θtrue))


surfind = surfaceindex(γ.I)
flush(stdout)
# The first guess for the tracer concentration should be close to the actual tracer concentration
# take first guess as θtrue
cvec=(vec(θtrue))
#θguess = unvec(θtrue,cvec)

ucnew = load("ucnew.jld")["ucnew"]
#first guess tracer control vector is near zero, and we want this to remain relatively small
u1 =  Field(0.0.*ones(size(γ.wet)),γ,θtrue.name,θtrue.longname,θtrue.units)
uvec = ucnew
u =  unvec(u1,ucnew)
q = zeros(size(uvec))
q[surfind] = cvec[surfind]

#test:set deep ocean to be constant temp
#cvec = ones(size(uvec))*0.5
#cvec[surfind] = q[surfind]


#We need an error covariance matrix
W⁻ = Diagonal(1 ./( ones(sum(γ.wet))).^2)#(1/sum(γ.wet))


Qdiag = 1 ./ ones(sum(γ.wet))
Q⁻ = Diagonal(Qdiag)

ind_1 = load("ind_1.jld")["ind_1"]
ind_2 = load("ind_2.jld")["ind_2"]
ufnew = load("ufnew.jld")["ufnew"]


Adummy = sparse(ind_1, ind_2, ones(length(ind_1))*0.1)


#test: apply A0 a bunch of times
#c_oneloc = zeros(length(cvec))
#c_oneloc[4000] = 15
#c_post = A0 * (A0 * (A0 * (A0 * (A0 * (A0 * (A0 * (A0 * c_oneloc)))))))
#c1_3D = unvec(θguess,Float32.(c_oneloc))
#c2_3D = unvec(θguess,Float32.(c_post))


# plot the difference
#level = 3 # your choice 1-33
#depth = γ.depth[level]

#cntrs = 0:0.5:15
#label = "Initial θ"
#planviewplot(c1_3D, depth, cntrs, titlelabel=label,cenlon=40)
#readline()


#cntrs = 0:0.5:15
#label = "Final θ"
#planviewplot(c2_3D, depth, cntrs, titlelabel=label,cenlon=40)
#readline()


non_zero_indices1, non_zero_indices2, non_zero_values = findnz(Adummy)

Afguess = sparse(non_zero_indices1, non_zero_indices2, ufnew)

non_zero_indices1, non_zero_indices2, non_zero_values = findnz(Afguess)

non_zero_indices = hcat(non_zero_indices1, non_zero_indices2)


convec = [uvec; non_zero_values]
ulength=length(uvec)

print(ulength)
print(sum(γ.wet))

# get sample J value
F = costfunction_gridded_model(convec,non_zero_indices,u,Adummy,ctrue,cvec,q,W⁻,Q⁻,γ)
fg!(F,G,x) = costfunction_gridded_model!(F,G,x,non_zero_indices,u,Adummy,ctrue,cvec,q,W⁻,Q⁻,γ)
fg(x) = costfunction_gridded_model(x,non_zero_indices,u,Adummy,ctrue,cvec,q,W⁻,Q⁻,γ)
f(x) = fg(x)[1]
J₀,gJ₀ = fg(convec)

#### gradient check ###################
# check with forward differences
ϵ = 1e-3
println(size(length(convec)))
#ii = rand(1:length(convec))
#ii=748529+10000
ii=4000
println("Location for test =",ii)
flush(stdout)
δu = copy(convec); δu[ii] += ϵ
∇f_finite = (f(δu) - f(convec))/ϵ
println(∇f_finite)

fg!(J₀,gJ₀,(convec+δu)./2) # J̃₀ is not overwritten
∇f = gJ₀[ii]
println(∇f)

# error less than 10 percent?
println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
flush(stdout)
#### end gradient check #################

print(length(convec))
# filter the data with an Optim.jl method
iterations = 20000
out = steadyclimatology_optim(convec,fg!,iterations)

# reconstruct by hand to double-check.
#ũ = unvec(θguess,Float32.(out.minimizer[begin:ulength]))

#println(typeof(ũ))

# reconstruct tracer map
#c₀ = θguess
#c̃  = θguess+ũ

#Δc̃ = c̃ - θtrue
#Δc₀ = θguess - θtrue

Anew = Adummy + sparse(non_zero_indices[:, 1], non_zero_indices[:, 2], out.minimizer[ulength+1:end])

maxA = maximum(Anew)
minA = minimum(Anew)

println("Max A:$maxA")
println("Min A:$minA")

ucout = out.minimizer[begin:ulength]

maxuc = maximum(ucout)
minuc = minimum(ucout)

println("Max uc:$maxuc")
println("Min uc:$minuc")


flush(stdout)
onesvec = ones(size(q))

#oldf = sum((non_zero_values).^2)
#newf = sum((out.minimizer[ulength+1:end]).^2)
tracer_cons1 = sum(((Adummy+Afguess)*(cvec+u) - q).^2)
tracer_cons2 = sum((Anew*(cvec+out.minimizer[begin:ulength]) - q).^2)
mass_cons1 = sum(((Adummy+ Afguess)*onesvec).^2)
mass_cons2 = sum((Anew*onesvec).^2)

println("old tracer cons:$tracer_cons1")
println("new tracer cons:$tracer_cons2")
println("old mass cons:$mass_cons1")
println("new mass cons:$mass_cons2")

b = getsurfaceboundary(θtrue)

save("ufnew.jld","ufnew", out.minimizer[ulength+1:end])
save("ucnew.jld","ucnew", out.minimizer[begin:ulength])





#A0


#Anew
#Alu2 = lu(Anew)


final_θ = steadyinversion(Anew,b,γ)
#

c1 = (vec(final_θ))

test = sum(Anew * c1 - q)


max_θ = maximum(filter(!isnan,final_θ.tracer))
min_θ = minimum(filter(!isnan,final_θ.tracer))

println("final_θ max:$max_θ")
println("final_θ min:$min_θ")
println("test:$test")

#writefield("final_theta.nc",final_θ)

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

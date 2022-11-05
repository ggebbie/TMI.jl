#=%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example : Find the distribution of a tracer given:              %
%       (a) the pathways described by A,                          %
%       (b) interior sources and sinks given by dC,               % 
%           that best fits observations, Cobs,                    %
%   and (c) inequality constraints on the tracer concentration.   %
%                                                                  %
% Mathematically, minimize J = (C-Cobs)^T W (C-Cobs) subject to    %
%                         AC = d + Gamma u                         %
%  where u is the estimated change in surface concentration.    % 
%
% See Supplementary Section 2, Gebbie & Huybers 2011.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% =#

import Pkg; Pkg.activate(".")

using Revise
using TMI
using Test
using GGplot
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# first guess of change to surface boundary conditions
# ocean values are 0
#u = zerosurfaceboundary(γ)
u = (;surface = zerosurfaceboundary(γ))
uvec = vec(u)
   
# take synthetic, noisy observations
y, W⁻, ctrue = synthetic_observations(TMIversion,"θ",γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
#b = getsurfaceboundary(y)
b = (;surface = getsurfaceboundary(y))

# get sample J value
F,G = costfunction_gridded_obs(uvec,Alu,b,u,y,W⁻,γ)
fg!(F,G,x) = costfunction_gridded_obs!(F,G,x,Alu,b,u,y,W⁻,γ)

fg(x) = costfunction_gridded_obs(x,Alu,b,u,y,W⁻,γ)
f(x) = fg(x)[1]
J₀,gJ₀ = fg(uvec)

#### gradient check ###################
# check with forward differences
ϵ = 1e-3
ii = rand(1:sum(γ.wet[:,:,1]))
println("Location for test =",ii)
δu = copy(uvec); δu[ii] += ϵ
∇f_finite = (f(δu) - f(uvec))/ϵ
println(∇f_finite)

fg!(J₀,gJ₀,(uvec+δu)./2) # J̃₀ is not overwritten
∇f = gJ₀[ii]
println(∇f)

# error less than 10 percent?
println("Percent error ",100*abs(∇f - ∇f_finite)/abs(∇f + ∇f_finite))
#### end gradient check #################

# filter the data with an Optim.jl method
iterations = 15
out = steadyclimatology(uvec,fg!,iterations)

# reconstruct by hand to double-check.
ũ = unvec(u,out.minimizer)

# apply to the boundary conditions
b̃ = adjustboundarycondition(b,ũ)

# reconstruct tracer map
c₀ = steadyinversion(Alu,b,γ)
c̃  = steadyinversion(Alu,b̃,γ)

Δc̃ = c̃ - ctrue
Δc₀ = c₀ - ctrue

# plot the difference
level = 15 # your choice 1-33
depth = γ.depth[level]
cntrs = -1:0.05:1
label = "Optimized misfit: Δc̃"
# Help: needs work with continents and labels
planviewplot(Δc̃, depth, cntrs, titlelabel=label) 

label = "First guess misfit: Δc₀"
planviewplot(Δc₀, depth, cntrs, titlelabel=label) 

# future step: box minimization to eliminate unreasonable temperatures.
# should be made more generic for other tracers
# From MATLAB:
#lbT = -2.*ones(Nsfc,1); % temperature lower bound: can not freeze.
#ubT = 40.*ones(Nsfc,1); % ad-hoc temperature upper bound: 40 C.

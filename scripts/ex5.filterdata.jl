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
# using Random, Statistics, LinearAlgebra, Dates #Shipped with Julia
# using Distributions, StatsBase #Core statistics
# using CSV, DataFrames #Basic Data
# using Plots, StatsPlots, LaTeXStrings, Measures #Plotting and Output
# using HypothesisTests, KernelDensity, GLM, Lasso, Clustering, Multivaria
# using Flux, Metalhead #Deep learning
# using Combinatorics, SpecialFunctions, Roots #Mathematical misc.
# using RDatasets, MLDatasets #Example datasets
# #uncomment if using R: using RCall #Inter

using Revise, TMI, GoogleDrive
using PyPlot, PyCall, Test
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "TMI_2010_2012_4x4x33"
A, Alu, γ, inputfile = config(TMIversion) 

# first guess of change to surface boundary conditions
# ocean values are 0
u₀ = zeros(Float64,sum(γ.wet[:,:,1]))

# take synthetic, noisy observations
y, W⁻, ctrue = sample_observations(TMIversion,"θ")

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
d₀ = tracerinit(γ.wet)
d₀[:,:,1] = y[:,:,1]

# check gradients in misfit_gridded_data!
fg(x) = misfit_gridded_data(x,Alu,y,d₀,W⁻,γ.wet)
J̃₀,gJ₀ = fg(u₀)
fg!(F,G,x) = misfit_gridded_data!(F,G,x,Alu,y,d₀,W⁻,γ.wet)

# check with forward differences
ϵ = 1e-5
ii = rand(1:sum(γ.wet[:,:,1]))
δu = copy(u₀); δu[ii] += ϵ
∇f_finite = (fg(δu)[1] - fg(u₀)[1])/ϵ # `[1]` to pick out cost function

fg!(J̃₀,gJ₀,(u₀+δu)./2) # J̃₀ is not overwritten
∇f = gJ₀[ii]

# error less than 10 percent?
@test (∇f - ∇f_finite)/abs(∇f + ∇f_finite) < 0.1

# sophisticated method: doesn't work due to issue #17
#gmisfitForward = x -> ForwardDiff.gradient(misfit, x); # g = ∇f
#gJforward = gmisfitForward(u₀)

# gradient check with autodiff? doesn't work: issue #17
#test = gradient(misfit,u₀)

# filter the data with an Optim.jl method
out = filterdata(u₀,Alu,y,d₀,W⁻,γ.wet)

# was cost function decreased?
@test out.minimum < J̃₀

# reconstruct by hand to double-check.
ũ = out.minimizer
J̃,gJ̃ = fg(ũ)
@test J̃ < J̃₀

# reconstruct tracer map
c₀ = steady_inversion(u₀,Alu,d₀,γ.wet)
c̃ = steady_inversion(ũ,Alu,d₀,γ.wet)

# view the surface
cntrs = -1:0.05:1

# what model depth level?
level = 15
Δc̃ = c̃[:,:,level] .- ctrue[:,:,level]
Δc₀ = c₀[:,:,level] .- ctrue[:,:,level]

figure()
contourf(γ.lon,γ.lat,Δc̃',cntrs)
#contour(γ.lon,γ.lat,Δc̃',cntrs)

figure()
contourf(γ.lon,γ.lat,Δc₀',cntrs)

# sample tests
#maximum(filter(!isnan,Δc₀))
#maximum(filter(!isnan,Δc̃))

# future step: box minimization to eliminate unreasonable temperatures.
# should be made more generic for other tracers
# From MATLAB:
#lbT = -2.*ones(Nsfc,1); % temperature lower bound: can not freeze.
#ubT = 40.*ones(Nsfc,1); % ad-hoc temperature upper bound: 40 C.

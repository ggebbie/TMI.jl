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

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# first guess of change to surface boundary conditions
# ocean values are 0
u₀ = zeros(Float64,sum(γ.wet[:,:,1]))

# take synthetic, noisy observations
y, W⁻, ctrue = synthetic_observations(TMIversion,"θ",γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
d₀ = tracerinit(γ.wet)
d₀[:,:,1] = y[:,:,1]

fg!(F,G,x) = costfunction_obs!(F,G,x,Alu,d₀,y,W⁻,γ)

# filter the data with an Optim.jl method
out = steadyclimatology(u₀,Alu,y,d₀,W⁻,fg!,γ)

# reconstruct by hand to double-check.
ũ = out.minimizer

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

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

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
d₀ = tracerinit(γ.wet)

# get observations at surface
# set them as surface boundary condition
PO₄obs = readtracer(TMIfile,"PO₄")
d₀[:,:,1] = PO₄obs[:,:,1]

# reconstruct tracer map
PO₄pre = steady_inversion(u₀,Alu,d₀,γ.wet)
PO₄ᴿ = regeneratedphosphate(TMIversion,Alu,γ)
PO₄total = PO₄pre + PO₄ᴿ

# invert for total phosphate in one step
qPO₄ = readtracer(TMIfile,"qPO₄")
dPO₄ = d₀ - qPO₄
PO₄direct = steady_inversion(u₀,Alu,dPO₄,γ.wet)

# view the surface
cntrs = 0:0.05:3.5

# what model depth level?
level = 15

figure()
contourf(γ.lon,γ.lat,PO₄recon[:,:,level]',cntrs)
contourf(γ.lon,γ.lat,PO₄total[:,:,level]',cntrs)
#contour(γ.lon,γ.lat,Δc̃',cntrs)

figure()
contourf(γ.lon,γ.lat,Δc₀',cntrs)

lon_section = 330; # only works if exact
Psection = section(PO₄total,lon_section,γ)
lims = 0:0.05:3.0
figure()
dyeplot(γ.lat,-γ.depth[33:-1:1],Psection[:,33:-1:1]', lims)


# oxygen distribution
O₂obs = readtracer(TMIfile,"O₂")
d₀[:,:,1] = O₂obs[:,:,1]
qO₂ = - 170 .* qPO₄
dO₂ = d₀ - qO₂
O₂recon = steady_inversion(u₀,Alu,dO₂,γ.wet)

# view the surface
cntrs = 0:10:400 # μmol/kg

# what model depth level?
level = 27

figure()
cf = contourf(γ.lon,γ.lat,O₂recon[:,:,level]',cntrs)
colorbar(cf)

lon_section = 202; # only works if exact
Psection = section(O₂recon,lon_section,γ)
figure()
dyeplot(γ.lat,-γ.depth[33:-1:1],Psection[:,33:-1:1]', cntrs)

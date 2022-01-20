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

using Revise, TMI, GoogleDrive
using PyPlot, PyCall, Test
#, Distributions, LinearAlgebra,  Zygote, ForwardDiff, Optim

TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
d = tracerinit(γ.wet)

δ¹³Cobs = readtracer(TMIfile,"δ¹³C")
                  
δ¹³Csfc = δ¹³Cobs[:,:,1] # nutrient surface boundary condition

# reconstruct tracer map
#δ¹³Cmodel = steady_inversion(Alu,d,γ.wet)

δ¹³Cmodel = steady_inversion(Alu,δ¹³Csfc,qPO₄,γ.wet)

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

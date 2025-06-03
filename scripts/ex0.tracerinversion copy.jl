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
cd("scripts")
import Pkg; Pkg.activate(".")

using Revise
using TMI
using GeoPythonPlot # will load optional extension
using Printf
using Distances
using LinearAlgebra
using BlockDiagonals

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

ones_like(x) = zero(x) .+ 1
lats = γ.lat .* ones_like(γ.lon)'
lons = ones_like(γ.lat) .* γ.lon' 
lons[(!).(γ.wet[:, :, 1])[:]] .= NaN
lats[(!).(γ.wet[:, :, 1])[:]] .= NaN

distances = distancematrix(γ)
@. gaussian_correlation(d) = exp(-0.5 * (abs(d)^2) / (500)^2) #decorrelation of 200km 

S_c = gaussian_correlation(distances)
S_c[S_c .== 1.0] .= 0.0
S_c[S_c .≈ 0.0] .= NaN
S_c = S_c .* 3

nx = 500
fig, ax = GeoPythonPlot.subplots()
pos_id = 1:length(lats[:])
cm = ax.pcolormesh(pos_id[1:nx], pos_id[1:nx], S_c[1:nx, 1:nx], rasterized = true, cmap = "Spectral_r")
fig.colorbar(cm)
fig

S_m = zero(lats[:]) .+ (0.3^2)
S_m = diagm(S_m)
S = Matrix(BlockDiagonal([S_c, S_m]))
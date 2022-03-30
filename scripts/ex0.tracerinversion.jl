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

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# get observations at surface
# set them as surface boundary condition
PO₄obs = readtracer(TMIfile,"PO₄")

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
bobs = getsurfaceboundary(PO₄obs,γ)

## preformed phosphate
PO₄pre = steadyinversion(Alu,bobs,γ)

## regenerated phosphate

qPO₄ = readtracer(TMIfile,"qPO₄")
q = getsource(qPO₄,γ)

# zero boundary condition
b₀ = zerosurfaceboundary(γ)
PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=q)
PO₄total = PO₄pre + PO₄ᴿ

## compute total phosphate directly
PO₄direct = steadyinversion(Alu,bobs,γ,q=q)

## how big is the maximum difference?
maximum(PO₄direct - PO₄total)
minimum(PO₄direct - PO₄total)

maximum(PO₄obs[γ.wet])

# view the surface
cntrs = 0:0.05:3.5

# what model depth level?
level = 15

figure()
contourf(γ.lon,γ.lat,PO₄ᴿ[:,:,level]',cntrs)
contourf(γ.lon,γ.lat,PO₄total[:,:,level]',cntrs)

lon_section = 330; # only works if exact
Psection = section(PO₄total,lon_section,γ)
lims = 0:0.05:3.0
figure()
sectionplot(γ.lat,γ.depth,Psection, lims)

# oxygen distribution
O₂obs = readtracer(TMIfile,"O₂")
d₀[:,:,1] = O₂obs[:,:,1]
qO₂ = - 170 .* qPO₄
dO₂ = d₀ - qO₂
O₂ᴿ = steady_inversion(u₀,Alu,dO₂,γ.wet)

# view the surface
cntrs = 0:10:400 # μmol/kg

# what model depth level?
level = 27

figure()
cf = contourf(γ.lon,γ.lat,O₂ᴿ[:,:,level]',cntrs)
colorbar(cf)

lon_section = 202; # only works if exact
Psection = section(O₂ᴿ,lon_section,γ)
figure()
#dyeplot(γ.lat,γ.depth[33:-1:1],Psection[:,33:-1:1], cntrs)
sectionplot(γ.lat,γ.depth,Psection, cntrs)

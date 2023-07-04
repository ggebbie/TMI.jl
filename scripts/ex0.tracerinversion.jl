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

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

# get observations at surface
# set them as surface boundary condition
yPO₄ = readfield(TMIfile,"PO₄",γ)

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.

# choice: BoundaryCondition or NamedTuple(BoundaryCondition)
#bPO₄ = getsurfaceboundary(yPO₄)
bPO₄ = (;surface = getsurfaceboundary(yPO₄))

## preformed phosphate
PO₄pre = steadyinversion(Alu,bPO₄,γ)

## read phosphate source
qPO₄ = readfield(TMIfile,"qPO₄",γ)

# zero boundary condition, choose one line of next two
#b₀ = zerosurfaceboundary(γ)
b₀ = (;surface = zerosurfaceboundary(γ,:PO₄,"phosphate","μmol/kg"))
PO₄ᴿ = steadyinversion(Alu,b₀,γ,q=qPO₄)
PO₄total = PO₄ᴿ + PO₄pre

## compute total phosphate directly
PO₄direct = steadyinversion(Alu,bPO₄,γ,q=qPO₄)

## how big is the maximum difference?
# could replace with abs function
@test maximum(PO₄direct - PO₄total) < 0.1
@test minimum(PO₄direct - PO₄total) > -0.1

## compare to observations
@test maximum(PO₄direct - yPO₄) < 1.2 # lenient test
@test minimum(PO₄direct - yPO₄) > -1.2

## Plot a plan view
# view the surface
cntrs = 0:0.2:4

# what model depth level?
level = 15
depth = γ.depth[level]
label = PO₄total.longname*", depth = "*string(depth)*" m"

# Help: needs work with continents and labels
planviewplot(PO₄total, depth, cntrs, titlelabel=label) 

## Plot a lat-depth section
lon_section = 330; # only works if exact
lims = 0:0.05:3.0
sectionplot(PO₄total,lon_section,lims)

## oxygen distribution
yO₂ = readfield(TMIfile,"O₂",γ)
bO₂ = getsurfaceboundary(yO₂)
# set the optional stoichiometric ratio.
O₂ = steadyinversion(Alu,bO₂,γ,q=qPO₄,r=-170.0)

# Plan view of oxygen at same depth as phosphate
# but different contours
cntrs = 0:10:400 # μmol/kg
planviewplot(O₂, depth, cntrs, titlelabel=label) 

# Section view of oxygen on same phosphate section.
sectionplot(O₂,lon_section,cntrs)

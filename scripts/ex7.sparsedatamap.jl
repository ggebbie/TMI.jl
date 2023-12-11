#=
     Find the distribution of a tracer given:
     (a) the pathways described by A or its LU decomposition Alu,
     (b) first-guess boundary conditions and interior sources given by d₀,
     (c) perturbations to the surface boundary condition u
    that best fits observations, y,
    according to the cost function,
    J = (ỹ - Ec)ᵀ W⁻¹ (ỹ - Ec)
    subject to Ac = d₀ + Γ u₀.                 
    W⁻¹ is a (sparse) weighting matrix.
    See Supplementary Section 2, Gebbie & Huybers 2011.
# Arguments
- `u`: control vector of surface tracer perturbations
- `Alu`: LU decomposition of water-mass matrix A
- `y`: observations on 3D grid
- `d₀`: first guess of boundary conditions and interior sources
- `W⁻`: weighting matrix best chosen as inverse error covariance matrix
- `wet`: BitArray mask of ocean points
=#
#using Revise, TMI, Interpolations, Statistics

import Pkg; Pkg.activate(".")

using Revise
using TMI
using Test
using GeoPythonPlot
using Interpolations
using Statistics
using LinearAlgebra

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

# first guess of change to surface boundary conditions
# how many randomly sampled observations?
N = 20

# first guess of change to surface boundary conditions
# ocean values are 0
#u = zerosurfaceboundary(γ)
u = (;surface = zerosurfaceboundary(γ))
uvec = vec(u)

# take synthetic, noisy observations
y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N)

# make a silly first guess for surface
#b = mean(y) * onesurfaceboundary(γ)
using Statistics
b = (;surface = mean(y) * onesurfaceboundary(γ))

# assume temperature known ± 5°C
σb = 5.0
Dg = gaussiandistancematrix(γ,σb,1000.0)
Q⁻ = inv(cholesky(Dg))

iters =5
                
out, f, fg, fg! = TMI.sparsedatamap(Alu,b,u,y,W⁻,wis,locs,Q⁻,γ,iters)

# reconstruct by hand to double-check.
ũ = unvec(u,out.minimizer)
b̃ = adjustboundarycondition(b,ũ)

# reconstruct tracer map
c₀ = steadyinversion(Alu,b,γ)
c̃  = steadyinversion(Alu,b̃,γ)

Δc̃ = c̃ - ctrue
Δc₀ = c₀ - ctrue

# plot the difference
level = 15 # your choice 1-33
depth = γ.depth[level]
cntrs = -10:0.5:10
label = "Optimized misfit: Δc̃"
# Help: needs work with continents and labels
planviewplot(Δc̃, depth, cntrs, titlelabel=label) 

label = "First guess misfit: Δc₀"
planviewplot(Δc₀, depth, cntrs, titlelabel=label) 

# % what is the uncertainty in the surface boundary condition?

# % how much did the data points reduce the error (globally).
# sqrt(sum((c-Tobs).^2)/Nfield)
# sqrt(sum((c0-Tobs).^2)/Nfield)

# % how much did the data points reduce the error (at data points).
# Nobs = length(y);
# sqrt(sum((E*c-y).^2)/Nobs)
# sqrt(sum((E*c0-y).^2)/Nobs)

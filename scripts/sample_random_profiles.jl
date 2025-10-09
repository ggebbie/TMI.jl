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
# using GeoPythonPlot
using Plots
using Interpolations
using Statistics
using LinearAlgebra

# Function to find the last non-NaN value in a vector
function last_non_nan(v)
    idx = findlast(!isnan, v)
    return idx !== nothing ? v[idx] : NaN
end

last_non_nan(m::Matrix) = last_non_nan.(eachcol(m))

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = TMI.config(TMIversion);

# first guess of change to surface boundary conditions
# how many randomly sampled observations?
N = 100

#Try with a Field
variable = "θ"
y, ctrue, ytrue, locs = random_profiles(TMIversion,variable,γ,N)

y_surface = y[1, :]
y_bottom = last_non_nan(y)

#Try with an array (in preparations for working with transient runs like GH19)
θtrue = readfield(TMIfile,variable,γ)
y, ctrue, ytrue, locs = random_profiles(θtrue.tracer,γ,N)
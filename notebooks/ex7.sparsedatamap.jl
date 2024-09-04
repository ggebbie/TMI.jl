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

import Pkg; Pkg.activate("scripts")

using Revise
using TMI
using Test
using GeoPythonPlot
import Pkg; Pkg.build("GeoPythonPlot")

using Interpolations
using Statistics
using LinearAlgebra

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion);

# first guess of change to surface boundary conditions
# how many randomly sampled observations?
N = 50

# first guess of change to surface boundary conditions
# ocean values are 0


# take synthetic, noisy observations
y, W⁻, ctrue, ytrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N)

p_list = [-0.05, 0.00, 0.05]
string(p_list[1])
c̃̄_dict = Dict()
for (i, p) in enumerate(p_list)
    nelm = (p == 0) ? 1 : N
    c̃̄_dict[p] = zeros(nelm)

    u = (;surface = zerosurfaceboundary(γ))
    uvec = vec(u)
    for j in 1:nelm
        yp = 1 .* y 
        yp[j] += p
        # make a silly first guess for surface
        b = (;surface = mean(yp) * onesurfaceboundary(γ))

        # assume temperature known ± 5°C
        σb = 5.0
        Dg = gaussiandistancematrix(γ,σb,1000.0)
        Q⁻ = inv(cholesky(Dg))
                        
        out, f, fg, fg! = TMI.sparsedatamap(Alu,b,u,yp,W⁻,wis,locs,Q⁻,γ)

        # reconstruct by hand to double-check.
        ũ = unvec(u,out.minimizer)
        b̃ = adjustboundarycondition(b,ũ)

        c̃  = steadyinversion(Alu,b̃,γ)
        # reconstruct tracer map
        c̃̄_dict[p][j]  = mean(c̃)
    end
end

W = (c̃̄_dict[p_list[end]] .- c̃̄_dict[p_list[1]]) ./ (p_list[end] - p_list[1])

θ =  readfield(TMIfile, "θ", γ)
mean(θ)

W' * y

c̃̄_dict[0.0]

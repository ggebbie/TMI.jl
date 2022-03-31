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
using Revise, TMI, Interpolations, Statistics

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# first guess of change to surface boundary conditions
# how many randomly sampled observations?
N = 20

# first guess of change to surface boundary conditions
# ocean values are 0
u = TMI.zerosurfaceboundary(γ)
uvec = u.tracer[u.wet]

# take synthetic, noisy observations
y, W⁻, ctrue, locs, wis = synthetic_observations(TMIversion,"θ",γ,N)

# make a silly first guess for surface
b = mean(y) * onesurfaceboundary(γ)

# assume temperature known ± 5°C
σb = 5.0
Q⁻ = 1.0/(σb^2)

# optimize the sparse data map with an Optim.jl method
out = sparsedatamap(u₀,Alu,d₀,y,W⁻,wis,locs,Q⁻,γ)

# missing: plots of results

# % what is the uncertainty in the surface boundary condition?
# err_d = 1;
# diagonal = false; % or set to false.
# if diagonal 
#   S     = 1./err_d.^2; 
# else
#   % if diagonal==false
#   % impose spatial smoothing in surface b.c.
#   lengthscale = 4; % horizontal lengthscale in units of gridcells
#   factor = 0.126.*(1/lengthscale)^2 ; 
#   load Del2_4deg.mat
#   S = factor.* sparse(1:Nsfc,1:Nsfc,1./err_d.^2);
#   S = S + Del2'*(lengthscale^4.*S*Del2);
# end
  

# % how much did the data points reduce the error (globally).
# sqrt(sum((c-Tobs).^2)/Nfield)
# sqrt(sum((c0-Tobs).^2)/Nfield)

# % how much did the data points reduce the error (at data points).
# Nobs = length(y);
# sqrt(sum((E*c-y).^2)/Nobs)
# sqrt(sum((E*c0-y).^2)/Nobs)

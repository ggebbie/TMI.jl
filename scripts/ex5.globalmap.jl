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
using Distributions, PyPlot, PyCall,
    LinearAlgebra, Zygote, ForwardDiff, Optim

url = "https://docs.google.com/uc?export=download&id=1Zycnx6_nifRrJo8XWMdlCFv4ODBpi-i7"
inputdir = "../data"

A, Alu, γ, inputfile = config(url,inputdir)

# take synthetic observations
# get observational uncertainty
# could be better to read indices to form wet mask rather than a sample variable.
θtrue = readtracer(inputfile,"θ")

#ΔPO₄ = readtracer(inputfile,"qpo4")
σθ = readtracer(inputfile,"σθ")

ntrue = tracerinit(γ.wet)
ntrue[γ.wet] = rand(Normal(),length(σθ[γ.wet])) .* σθ[γ.wet]

y = θtrue .+ ntrue

# get cost function (J) based on model misfit
# here the data-model misfit is weighted by the expected error

# weighting matrix
W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
Wⁱ = (1/sum(γ.wet)) .* Diagonal(1 ./σθ[γ.wet].^2)

# any penalty for varying controls?
#W = sum(γ.wet) .* Diagonal(σθ[γ.wet].^2)
Qⁱ = 0.0 # already penalized by surface data

# a first guess: observed surface boundary conditions are perfect.
# set surface boundary condition to the observations.
# below surface = 0 % no internal sinks or sources.
d₀ = tracerinit(γ.wet)
d₀[:,:,1] = y[:,:,1]

# first-guess reconstruction of observations
ỹ₀ = tracerinit(γ.wet)
ỹ₀[γ.wet] = Alu\d₀[γ.wet]

ñ₀ = tracerinit(γ.wet)
ñ₀ = y .- ỹ₀

# starting value of J̃
J̃₀ = ñ₀[γ.wet]'* Wⁱ * ñ₀[γ.wet]

# first guess of change to surface boundary conditions
# ocean values are 0
u₀ = Vector{eltype(u₀)}(undef,sum(γ.wet[:,:,1]))
fill!(u₀,zero(eltype(u₀vec)))

# can you make this a function of u?
J̃₀b,gJb = misfit_gridded_data(u₀,Alu,y,d₀,Wⁱ,Qⁱ,γ.wet)
isapprox(J̃₀,J̃₀b)

fg(x) = misfit_gridded_data(x,Alu,y,d₀,Wⁱ,Qⁱ,γ.wet)
F0,G0 = fg(u₀)
fg!(F,G,x) = misfit_gridded_data!(F,G,x,Alu,y,d₀,Wⁱ,Qⁱ,γ.wet)
F1,G1 = fg!(F0,G0,u₀)

#gmisfit(x) = gmisfit_gridded_data(x,Alu,y,d₀,Wⁱ,Qⁱ,γ.wet)

# use explicitly supplied gradient function
#gJ₀ = gmisfit_gridded_data(u₀,Alu,y,d₀,Wⁱ,Qⁱ,γ.wet)
#gJ₀b = gmisfit_gridded_data(u₀vec,Alu,y,d₀,Wⁱ,Qⁱ,γ.wet)
isapprox(gJ₀[γ.wet[:,:,1]],gJ₀b)

# check with forward differences
# elementary method
ϵ = 1e-5
δu = copy(u₀); δu[20,20] += ϵ
g_elementary = (misfit(δu) - misfit(u₀))/ϵ

# error less than 10 percent?
(gJ₀[20,20]-g_elementary)/abs(gJ₀[20,20]+g_elementary) < 0.1
    
# sophisticated method: doesn't work due to issue #17
gmisfitForward = x -> ForwardDiff.gradient(misfit, x); # g = ∇f
gJforward = gmisfitForward(u₀)

# gradient check with autodiff? doesn't work: issue #17
#test = gradient(misfit,u₀)

# optimize with Optim.jl

# no gradient
#optimize(misfit,u₀,NelderMead())

# with gradient
G = gmisfit(u₀vec)
function g!(G, x)
    G = gmisfit(x)
    return G
end
out = optimize(misfit,g!,u₀vec,LBFGS())

#= minimize J̃ to approximately one by modifying the surface boundary condition within its uncertainty.
 mathematical form: Find uT such that:
 J = (Tmod-Tobs)'*inv(WT)*(Tmod-Tobs) is minimized 
subject to:  A*Tmod = dT + Gamma * uT. =#


# Γ is a function that takes surface values and returns a 3D field with zeroes elsewhere
#test =  Γ(u₀,γ.wet)

lbT = -2.*ones(Nsfc,1); % temperature lower bound: can not freeze.
ubT = 40.*ones(Nsfc,1); % ad-hoc temperature upper bound: 40 C.

%% 4 numerical methods: 1) Constrained minimization with Lagrange multipliers 
%(but without inequality constraints), 2) quadratic programming using full Hessian,
% 3) quadratic programming using functional form of Hessian, 
% 4) Unconstrained minimization according to Gebbie et al. 2015 (QSR). 

% Here we proceed with method #1. Warning: convergence may take 30+ minutes.
options = optimset('Algorithm','interior-point','Display','iter', ...
                  'GradObj','on','LargeScale','on','TolX',1e-1);
% also can try 'Algorithm','trust-region-reflective'
noncons = 0;
isfc = find(kt==1);
uTtilde= fmincon(@(x)objfun(x,A,invWT,Tobs,isfc,inotmixlyr,noncons),Tobs(isfc),[],[],[],[],lbT,ubT,[],options);
J2 = objfun(uTtilde,A,invWT,Tobs,isfc,inotmixlyr,noncons)./N % expect J2 ~ 1, J2<J1
d2 = zeros(N,1); d2(isfc) = uTtilde;
Tmod2 =  Q * (U \ (L \ (P * (R \ d2)))) ; % best estimate
=#

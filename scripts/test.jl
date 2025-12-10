import Pkg; Pkg.activate(".")

using Revise
using TMI
using LinearAlgebra

TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion)
θ = readfield(TMIfile, "θ", γ)

# Define the standard deviation of the temperature data
θ_error = 0.1

# Calculate the variance
θ_variance = θ_error^2

# Create a precision matrix (diagonal in this case)
# The precision is the inverse of the variance
precision = Diagonal((1/θ_variance) .* ones(length(vec(θ_obs))))

# Create the Observation object
c_θ = Observations(θ; W=precision)



obs_values = [2.5, 3.1, 1.8] # e.g., tracer concentrations

# 2. Define the locations of the observations [lon, lat, depth]
locs = [
    [-150.0, 30.0, 500.0],  # Location 1
    [-150.0, 35.0, 1000.0], # Location 2
    [-150.0, 40.0, 1500.0]  # Location 3
]

# 3. Specify the uncertainty of the observations
obs_error = 0.2
obs_variance = obs_error^2
precision = Diagonal((1/obs_variance) .* ones(length(obs_values)))

# 4. Create the Observations object
sparse_obs = Observations(obs_values;
                         locs=locs,
                         γ=γ,
                         W=precision)

# 5. Compare to model at observation sites
n = observe(θ, sparse_obs, γ) .-sparse_obs.values


τ = 1 / 2

sparse_obs_decay = Observations(obs_values;
                               locs=locs,
                               γ=γ,
                               W=precision,
                               decay_rate=τ)
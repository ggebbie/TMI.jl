# Observations

An `Observation` is a data structure that holds observed tracer data and its associated uncertainty. 
This is an important input for inverse modeling, where the goal is to find model parameters that best match the observations.

## Creating an `Observation`

An `Observation` is created using the `Observations` constructor, which takes the observed data and a precision matrix `W` as input. The precision matrix is the inverse of the covariance matrix of the data, and it is used to weight the cost function in the inversion.

Here's an example of how to create an `Observation` object for a temperature field:

```julia
# Load a gridded temperature data
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
```

In this example:

1.  We first load the observed temperature data `θ_obs`.
2.  We define the uncertainty of the data as a standard deviation `θ_error`.
3.  We calculate the variance `θ_variance`.
4.  We create a diagonal precision matrix. A diagonal precision matrix assumes that the errors in the data are uncorrelated.
5.  We then create the `Observation` object `c_θ` with the observed data and the precision matrix.

This `Observation` object can then be used in the inversion to find the model parameters that best fit the data.

## Sparse Observations

The `Observations` constructor can also handle sparse observations, where data is only available at a limited number of locations. This is done by providing a `locs` vector, which specifies the locations of the observations.

Here is an example of how to create an `Observations` object for sparse data and compare it to a gridded field:

```julia
# 1. Define the sparse observation values
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

```

In this case:
- `obs_values` is a vector of the tracer values.
- `locs` is a vector of coordinates, where each coordinate is a vector of `[longitude, latitude, depth]`.
- `γ` is the model grid, which is required for interpolating the model to the observation locations.
- `W` is the precision matrix, representing the uncertainty of the sparse observations.

## Sparse Observations with Decay Rate

A `decay_rate` can be provided for tracers that decay over time (e.g., radiocarbon).

Here is an example that builds on the sparse observation example by adding a decay rate:

```julia
# (assuming obs_values, locs, γ, and precision are defined as in the previous example)

# 1. Define the decay rate (e.g., for a radioactive tracer)
#    The decay rate should have units of 1/s.
τ = 1 / 2

# 2. Create the Observations object with a decay rate
sparse_obs_decay = Observations(obs_values;
                               locs=locs,
                               γ=γ,
                               W=precision,
                               decay_rate=τ)
```

In this case, the `Observations` constructor will use the `decay_rate` `τ`. 
This `sparse_obs_decay` object can then be used in an inversion to constrain the model with the sparse, decaying tracer data by modifying the watermassmatrix `A` via `A - sparse_obs_decay.decay_rate * I.` .

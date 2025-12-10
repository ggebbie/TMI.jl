# Controls

In the context of `TMI.jl`'s inverse modeling capabilities, "controls" refer to the parameters of the model that are adjusted to best fit the observations. These can include surface boundary conditions, interior sources/sinks of tracers, and water mass fractions.

The `TMI.setup_inversion` function is used to define these controls, along with their prior estimates, uncertainties (variances), and optimization bounds.

## Setting up Controls

The `TMI.setup_inversion` function takes keyword arguments for `boundary`, `source`, and `mass_fraction`. Each of these arguments is a `NamedTuple` that can contain:

-   `prior`: The initial estimate for the control variable.
-   `variance`: The variance of the prior estimate, representing our confidence in the prior.
-   `lower_bound`: The lower bound for the control variable during optimization.
-   `upper_bound`: The upper bound for the control variable during optimization.

Here's an example of how to set up controls for a simple inversion:

```julia
#load in gridded fields
TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

cobs = (θ =  readfield(TMIfile, "θ", γ), S = readfield(TMIfile, "Sp", γ))  

# Priors for surface boundary conditions are taken from observations
u₀ = map(v -> getsurfaceboundary(v), cobs)

# If you want to constrain a new boundary condition according to a 
# prior, you must specify a prior uncertainty 
u₀_uncertainty_variance = (θ = 0.1, S = 0.01)

# Priors for interior sources are set to nothing (no sources)
q₀ = (θ = nothing, S = nothing)

# Isotropic mass fractions as a prior estimate
m0 = massfractions_isotropic(γ)

# Define optimization bounds for the boundary control variables
u_lower = (θ = -2.0, S = 0.0)
u_upper = (θ = 35.0, S = 45.0)

# Set up the controls using setup_inversion
controls = TMI.setup_inversion(γ;
    boundary = (
        prior = u₀,
        variance = u₀_uncertainty_variance,
        lower_bound = u_lower,
        upper_bound = u_upper
    ),
    source = (
        prior = q₀
        # no variance or bounds needed if there are no sources
    ),
    mass_fraction = (
        prior = m0,
        variance = 1.0,
        lower_bound = map(v -> 0.0, m0),
        upper_bound = map(v -> 1.1, m0)
    )
)
```

In this example:

1.  We define priors for the boundary conditions (`u₀`), sources (`q₀`), and mass fractions (`m0`).
2.  We define lower and upper bounds for the boundary conditions.
3.  We then pass all this information to `TMI.setup_inversion` to create the `controls` object.

This `controls` object, along with the `Observations`, is then passed to the optimization routine to find the optimal set of control parameters.

# Controls

In the context of `TMI.jl`'s inverse modeling capabilities, "controls" refer to the parameters of the model that are adjusted to best fit the observations. These can include surface boundary conditions, interior sources/sinks of tracers, and water mass fractions.

Use the `BoundaryControls`, `SourceControls`, and `MassFracControls` constructors to define priors, uncertainties (variances or covariances), and optimization bounds for each component. Assemble them into a single `Controls` object with `TMI.Controls`, which handles vectorization and bounds for the optimizer.

## Setting up Controls

Each control constructor accepts optional keyword arguments for `variance` or `covariance`, `lower_bound`, `upper_bound`, and (for sources) `dependencies`. If you pass `nothing` or an empty `NamedTuple` for a component, the corresponding control block is skipped.

Here's an example of how to set up controls for a simple inversion using the current API:

```julia
using TMI, LinearAlgebra
TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

cobs = (θ =  readfield(TMIfile, "θ", γ), S = readfield(TMIfile, "Sp", γ))  

u₀ = map(v -> getsurfaceboundary(v), cobs)

tracer_error = (θ = 0.1, S = 0.01)
tracer_error_variance = map(x -> x^2, tracer_error)

q₀ = (θ = nothing, S = nothing) # no interior sources
m0 = massfractions_isotropic(γ)  # isotropic mass fractions

u_lower = (θ = -2.0, S = 0.0)
u_upper = (θ = 35.0, S = 45.0)

boundary_controls = BoundaryControls(
    u₀;
    variance = tracer_error_variance,
    lower_bound = u_lower,
    upper_bound = u_upper
);

source_controls = SourceControls(q₀); # no variance/bounds needed when everything is fixed

massfrac_controls = MassFracControls(
    m0;
    variance = 1.0,
    lower_bound = map(v -> 0.0, m0),
    upper_bound = map(v -> 1.1, m0),
    γ = γ
);

controls = TMI.Controls(γ;
    boundary = boundary_controls,
    source = source_controls,
    massfrac = massfrac_controls
);

#returns the controls vector. in this case, 
#these are the θ and S boundary conditions 
#as well as mass fractions, all to be optimized. 
vec(controls); 
```

In this example:

1.  Priors are defined for boundary conditions (`u₀`) and mass fractions (`m0`); sources are disabled by setting `q₀` fields to `nothing`.
2.  Variances and bounds are supplied directly to each constructor.
3.  The component controls are combined into a single `Controls` object, which provides vectorized parameters and bounds for your optimizer.

This `controls` object, along with `Observations`, is then passed to the optimization routine to find the optimal set of control parameters.

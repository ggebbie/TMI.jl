import Pkg; Pkg.activate(".")

using Revise
using LinearAlgebra
using TMI
using Statistics
using Test
using Optim, LineSearches, BenchmarkTools
using LinearSolve, IncompleteLU

# Helper function to calculate percent difference
percent_difference(x, y) = @. 100 * (x - y) / y

# Load a pre-configured TMI model and grid.
# This sets up the transport matrix `A`, its LU factorization `Alu`,
# the grid `γ`, and other model parameters.
TMIversion = "modern_180x90x33_GH11_GH12"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

# Load observed tracer fields from the TMI NetCDF file.
# These will be the target for the inversion.
cobs = (θ =  readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
    # Other tracers can be included here if available.
    # δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
    # P★ = preformedphosphate(TMIversion,Alu,γ),
    # δ¹³C★ = TMI.preformedcarbon13(TMIversion,Alu,γ)
)

# Use the loaded observations as the "true" state for this inversion experiment.
c0 = cobs

# Isotropic mass fractions as a prior estimate.
m0 = massfractions_isotropic(γ)

# Define tracer standard deviations (σ) for boundary conditions and observations.
# These represent our confidence in the prior values.
tracer_stddev = (θ = 0.01, S = 0.001)

# Square them to get variances (σ²), which are used in the cost function.
tracer_variance = map(x -> x^2, tracer_stddev)

# Priors for surface boundary conditions (u₀) are taken from the observations.
u₀ = map(v -> getsurfaceboundary(v), cobs)

# Priors for interior sources (q₀) are set to nothing (no sources).
q₀ = map(v -> nothing, cobs)

# Define optimization bounds for the control variables.
u_lower = (θ = -2.0, S = 0.0)
u_upper = (θ = 35.0, S = 45.0)

# This bundles all priors, uncertainties, and bounds into a single object
# for the optimization routine, using the new `setup_inversion` API.
controls = TMI.setup_inversion(γ;
    boundary = (
        prior = u₀, 
        variance = tracer_variance,
        lower_bound = u_lower,
        upper_bound = u_upper
    ),
    source = (prior = q₀,), # No variance specified for source in original
    mass_fraction = (prior = m0, variance = 1.0) # Unit variance for mass fractions
)

# Create observation objects with their corresponding precision (inverse variance).
# The precision matrix W is 1/σ², weighting the cost function.
observation_precision = map((v, σ²) -> Diagonal( (1/σ²) .* one.(vec(v)) ), c0, tracer_variance)
c_obs = map((v, w) -> Observations(v; W = w), c0, observation_precision)

"""
Sets up and runs the L-BFGS optimization for the TMI inversion.
"""
function constrained_global_optimization(controls::ControlParameters, c_obs, γ::Grid; 
                                        x0::Union{Nothing,Vector{T}} = nothing) where T

    # The `fg!` function calculates the objective (F) and gradient (G) for Optim.
    function fg!(F, G, x)
        optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=nothing)
    end

    # If no initial guess is provided, create a naive one from the priors.
    if isnothing(x0)
        x0 = vcat([vec(controls.u₀), vec(controls.q₀), randn(length(vec(controls.m₀)))]...)
    end

    # Run the L-BFGS optimization.
    result_opt_fg = Optim.optimize(
        Optim.only_fg!(fg!), x0,
        LBFGS(;m = 10, 
                alphaguess = LineSearches.InitialHagerZhang(α0=NaN), 
                linesearch = LineSearches.HagerZhang()),
        Optim.Options(f_abstol = 1e-12, g_tol = 1e-8,
                    iterations = 2,
                    store_trace = false,
                    show_trace = true, show_warnings = false));
    return result_opt_fg
end

# Construct the initial guess vector from the priors.
x0 = vcat([vec(controls.boundary.u₀), vec(controls.source.q₀), vec(controls.massfrac.m₀) ]...)

# Run the inversion.
result_opt_fg = constrained_global_optimization(controls, c_obs, γ; x0 = x0)

# Uncomment to benchmark the optimization.
# @btime constrained_global_optimization(controls, c_obs, γ; x0 = x0);
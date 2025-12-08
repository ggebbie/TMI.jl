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

# Define standard deviations (σ) for the boundary condition priors.
# These values are taken directly from the manual setup in total_inversion_3d.jl
w = (θ = 0.01, S = 0.001)

# This bundles all priors, uncertainties, and bounds into a single object
# for the optimization routine, using the `setup_inversion` API,
# but configured to match the manual setup of total_inversion_3d.jl.
controls = TMI.setup_inversion(γ;
    boundary = (
        prior = u₀,
        variance = w # NOTE: Passing stddev as variance to match total_inversion_3d.jl's Qᵤ of 1/σ
    ),
    source = (
        prior = q₀,
    ),
    mass_fraction = (
        prior = m0,
        variance = 1.0 # Matches Qₘ = Diagonal(one.(vec(m0)))
    )
)

# Create observation objects with unit weights, as in total_inversion_3d.jl
W_full = map(v -> Diagonal(one.(vec(v))), c0)
c_obs = map((v, w) -> Observations(v; W = w), c0, W_full)

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
        x0 = vcat([vec(controls.ub), vec(controls.uq), randn(length(vec(controls.m)))]...)
    end

    # Run the L-BFGS optimization.
    result_opt_fg = Optim.optimize(
        Optim.only_fg!(fg!), x0,
        LBFGS(;m = 10, 
                alphaguess = LineSearches.BackTracking(), 
                linesearch = LineSearches.HagerZhang()),
        Optim.Options(f_abstol = 1e-12, g_tol = 1e-8,
                    iterations = 15, outer_iterations = 2, 
                    store_trace = false,
                    show_trace = true, show_warnings = false));
    return result_opt_fg
end

# Construct the initial guess vector from the priors.
x0 = vcat([vec(controls.boundary.ub), vec(controls.source.uq), randn(length(vec(controls.massfrac.m)))]...)

function fg!(F, G, x)
    optim_fg_constrained_global_costfunction!(F, G, x, controls, c_obs, γ; locs=nothing)
end


# Run the inversion.
result_opt_fg = constrained_global_optimization(controls, c_obs, γ; x0 = x0)

# Uncomment to benchmark the optimization.
# @btime constrained_global_optimization(controls, c_obs, γ; x0 = x0)
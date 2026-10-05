#=
Invert for the water-mass fractions that fit six TMI tracers, holding the
surface boundary conditions and the phosphate source at their first guess.
Gradients come from Enzyme and the optimizer is Ipopt.

Run from the repository with Julia 1.13, with Enzyme, Optimization, and
OptimizationIpopt in the global environment:
    JULIA_LOAD_PATH="@:@v1.13:@stdlib" julia +1.13 --project=. scripts/mass_fraction_inversion.jl
Plot the result with scripts/inversion_diagnostics.jl.
=#
using Enzyme
using Optimization
using OptimizationIpopt
using TMI

TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

# observations: TMI tracers with uniform uncertainty
c = (θ = readfield(TMIfile, "θ", γ),
    S = readfield(TMIfile, "Sp", γ),
    δ¹⁸O = readfield(TMIfile, "δ¹⁸Ow", γ),
    PO₄ = readfield(TMIfile, "PO₄", γ),
    NO₃ = readfield(TMIfile, "NO₃", γ),
    O₂ = readfield(TMIfile, "O₂", γ))
σ = (θ = 0.1, S = 0.01, δ¹⁸O = 0.2, PO₄ = 0.05, NO₃ = 1.0, O₂ = 5.0)
y = Observations(c, σ, γ)

# first guess: observed surface values, TMI phosphate source, isotropic mass fractions
b₀ = map(getsurfaceboundary, c)
q₀ = (qPO₄ = readsource(TMIfile, "qPO₄", γ),)
m₀ = massfractions_isotropic(γ)

# stoichiometric ratios for A = watermassmatrix(m, γ), whose interior diagonal
# is +1: signs are opposite to ex0, which uses the A of the TMI file
r = (PO₄ = (qPO₄ = -1.0,), NO₃ = (qPO₄ = -15.5,), O₂ = (qPO₄ = 170.0,))

# adjust only the mass fractions, within [0, 1]
controls = (m = (σ = 1.0, lb = 0.0, ub = 1.0),)
inversion = Inversion(y; b₀, q₀, r, m₀, controls)

optimizer = IpoptOptimizer(
    hessian_approximation = "limited-memory",
    limited_memory_max_history = 27,
    acceptable_tol = 1.0e-6,
    mu_strategy = "adaptive",
    adaptive_mu_globalization = "kkt-error",
    nlp_scaling_method = "gradient-based",
)
runinversion(inversion, optimizer; name = "mass_fraction_inversion",
    iterations = 10, checkpoint_interval = 1)

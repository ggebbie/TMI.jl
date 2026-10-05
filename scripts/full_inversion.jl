#=
Invert for the surface boundary conditions, the phosphate source, and the
water-mass fractions that fit gridded WOCE observations of θ, S, δ¹⁸O, PO₄,
NO₃, and O₂. Gradients come from Enzyme and the optimizer is Ipopt.
The WOCE files are read from a sibling checkout of regridding_WOCE_for_TMI.

Run from the repository with Julia 1.13, with Enzyme, Optimization, and
OptimizationIpopt in the global environment:
    JULIA_LOAD_PATH="@:@v1.13:@stdlib" julia +1.13 --project=. scripts/full_inversion.jl
Plot the result with scripts/inversion_diagnostics.jl.
=#
using Enzyme
using Optimization
using OptimizationIpopt
using TMI

TMIversion = "modern_90x45x33_G14_v2"
A, Alu, γ, TMIfile, L, B = config(TMIversion);

WOCEdir = normpath(TMI.pkgdir("..", "regridding_WOCE_for_TMI", "data"))
WOCEfile = joinpath(WOCEdir, "TMI_gridded_WOCE_4x4_Variables_33_levels.nc")
WOCEerrorfile = joinpath(WOCEdir, "TMI_gridded_WOCE_Errors_4x4_Variables_33_levels.nc")

# WOCE variable and TMI name of each tracer
WOCEnames = (θ = "potential_temperature", S = "salinity", δ¹⁸O = "d18o",
    PO₄ = "phosphate", NO₃ = "nitrate", O₂ = "oxygen")
TMInames = (θ = :θ, S = :Sₚ, δ¹⁸O = :δ¹⁸Ow, PO₄ = :PO₄, NO₃ = :NO₃, O₂ = :O₂)

# observations, their standard deviations, and the horizontal decorrelation
# length [km] of their errors
c = map((variable, name) -> readfield(WOCEfile, variable, γ; name), WOCEnames, TMInames)
σ = map((variable, name) -> readfield(WOCEerrorfile, variable, γ; name = Symbol("σ", name)),
    WOCEnames, TMInames)
Ly = readfield(WOCEerrorfile, "decorrelation_length", γ;
    name = :L, longname = "decorrelation length", units = "km")
y = Observations(c, σ, γ; L = Ly)

# first guess: observed surface values, a uniform phosphate source, isotropic mass fractions
b₀ = map(getsurfaceboundary, c)
q₀ = (qPO₄ = 3.0e-4 * onesource(γ, :qPO₄, "local source of phosphate", "μmol/kg"),)
m₀ = massfractions_isotropic(γ)

# stoichiometric ratios for A = watermassmatrix(m, γ), whose interior diagonal
# is +1: signs are opposite to ex0, which uses the A of the TMI file
r = (PO₄ = (qPO₄ = -1.0,), NO₃ = (qPO₄ = -15.5,), O₂ = (qPO₄ = 170.0,))

# physical limits on each surface boundary condition
lb = (θ = -2.0, S = 0.0, δ¹⁸O = -10.0, PO₄ = 0.0, NO₃ = 0.0, O₂ = 0.0)
ub = (θ = 35.0, S = 45.0, δ¹⁸O = 10.0, PO₄ = 10.0, NO₃ = 45.0, O₂ = 500.0)

# adjustments: standard deviation σ, smoothing length L [km], and physical limits lb, ub.
# The phosphate source is adjusted on a log scale, qPO₄ = qPO₄₀ eᵘ, which keeps it
# positive; its σ is therefore in log units (log 50: one standard deviation is a
# factor of 50), while lb and ub stay limits on qPO₄ in its own units.
σb = map(getsurfaceboundary, σ)
controls = (b = map((σₖ, lbₖ, ubₖ) -> (σ = σₖ, L = 1_000.0, lb = lbₖ, ub = ubₖ), σb, lb, ub),
    q = (qPO₄ = (σ = log(50.0), lb = 6.0e-6, ub = 1.5e-2, logscale = true),),
    m = (σ = 1.0, lb = 0.0, ub = 1.0))
inversion = Inversion(y; b₀, q₀, r, m₀, controls)

optimizer = IpoptOptimizer(
    hessian_approximation = "limited-memory",
    limited_memory_max_history = 10,
    acceptable_tol = 1.0e-6,
    mu_strategy = "adaptive",
    adaptive_mu_globalization = "kkt-error",
    nlp_scaling_method = "gradient-based",
)
runinversion(inversion, optimizer; name = "full_inversion",
    iterations = 6000, checkpoint_interval = 500)
